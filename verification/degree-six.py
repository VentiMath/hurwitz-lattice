"""Exact verification of the degree-six computation in Appendix B.

Run with Python 3 and SymPy. Checks integral quaternionic action matrices,
the mod-two fixed space, every entry of the printed six-orbit certificate,
and the determinant of the resulting 24 orbit vectors. The mathematical
proof in the manuscript does not depend on this script.
"""
import argparse
from pathlib import Path
from itertools import combinations_with_replacement, product


def main():
    parser = argparse.ArgumentParser(
        description="Verify the degree-six computation in Appendix B exactly."
    )
    parser.add_argument(
        "--matrices", type=Path, metavar="PATH",
        help="write exact matrices to PATH after every check passes (overwrites PATH)",
    )
    args = parser.parse_args()
    if not __debug__:
        parser.error("verification requires assertions; do not run with python -O or -OO")

    import sympy as s

    # D4 basis in degree-two homology coordinates a,b,c,d.
    B4=s.Matrix([[2,-1,-1,-1],[0,1,0,0],[0,0,1,0],[0,0,0,1]])
    I=s.Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]])
    J=s.Matrix([[0,0,-1,0],[0,0,0,1],[1,0,0,0],[0,-1,0,0]])
    triples=list(combinations_with_replacement(range(4),3))
    def sym3(M):
        T=s.zeros(20)
        for col,abc in enumerate(triples):
            for rst in product(range(4),repeat=3):
                T[triples.index(tuple(sorted(rst))),col]+=s.prod(M[r,a] for r,a in zip(rst,abc))
        return T

    def rref2(M):
        A=[[int(x)%2 for x in row] for row in M.tolist()]
        piv=[];row=0
        for col in range(M.cols):
            found=next((r for r in range(row,len(A)) if A[r][col]),None)
            if found is None:continue
            A[row],A[found]=A[found],A[row]
            for r in range(len(A)):
                if r!=row and A[r][col]:A[r]=[x^y for x,y in zip(A[r],A[row])]
            piv.append(col);row+=1
        null=[]
        for f in range(M.cols):
            if f in piv:continue
            v=[0]*M.cols;v[f]=1
            for r,p in enumerate(piv):v[p]=A[r][f]
            null.append(v)
        return piv,s.Matrix.hstack(*(s.Matrix(v) for v in null)) if null else s.zeros(M.cols,0)

    a,b,c,d=s.symbols('a b c d'); ms=[a,b,c,d]
    C=s.zeros(20,4)
    corrections=[-4*a**3,-a*(b-a)**2,-a*(c-a)**2,-a*(d-a)**2]
    for j,poly in enumerate(corrections):
        p=s.Poly(poly,*ms)
        for row,abc in enumerate(triples):
            C[row,j]=p.coeff_monomial(s.prod(ms[t] for t in abc))
    B=s.zeros(24)
    B[:20,:20]=sym3(B4)
    B[:20,20:]=C
    B[20:,20:]=B4
    assert B.det()==2**16
    Ii=B.inv()*s.diag(sym3(I),I)*B
    Jj=B.inv()*s.diag(sym3(J),J)*B
    assert all(x.q==1 for x in Ii)
    assert all(x.q==1 for x in Jj)
    assert Ii**2==Jj**2==-s.eye(24)
    assert Ii*Jj==-Jj*Ii
    fixed=(Ii-s.eye(24)).col_join(Jj-s.eye(24))
    piv,ker=rref2(fixed)
    print('Triple order in (v1,t1^(2),t1^(3),t1^(4)):',triples)
    print('Hurewicz determinant:',B.det())
    print('i off-diagonal block (cubic part of i acting on the four degree-six generators):')
    print(Ii[:20,20:])
    print('j off-diagonal block:')
    print(Jj[:20,20:])
    print('Mod-two fixed-space dimension:',ker.cols)
    print('Fixed basis columns:',ker)
    print('Images of these fixed vectors in quotient D/2D:',ker[20:,:])
    print('Image dimension in quotient:',len(rref2(ker[20:,:])[0]))
    print('All action entries integral, all quaternion relations verified.')
    # Compact freeness certificate. In this paragraph a,b,c,d are formal coordinates
    # for (v1,t1^(2),t1^(3),t1^(4)), not homology coordinates.
    f,g,u,v=a,b+c+d,b,c
    p_polys=[f*u*v,g*u*v,u**3,u*u*v]
    def homotopy_cubic_column(poly):
        pol=s.Poly(poly,*ms)
        return s.Matrix([pol.coeff_monomial(s.prod(ms[t] for t in abc)) for abc in triples]+[0]*4)
    gens=[homotopy_cubic_column(p) for p in p_polys]+[s.eye(24)[:,21],s.eye(24)[:,22]]
    mod=lambda m:m.applyfunc(lambda x:int(x)%2)
    xy=mod((Ii-s.eye(24))*(Jj-s.eye(24)))
    F,G,U,V=s.symbols('f g u v')
    adapted_vars=(F,G,U,V)
    adapted_norms=[]
    for gen in gens:
        vec=mod(xy*gen)
        assert vec[20:,:]==s.zeros(4,1)
        poly=sum(vec[r]*s.prod(ms[t] for t in abc) for r,abc in enumerate(triples))
        pol=s.Poly(s.expand(poly.subs({a:F,b:U,c:V,d:G+U+V},simultaneous=True)),*adapted_vars)
        poly=sum((int(co)%2)*s.prod(t**e for t,e in zip(adapted_vars,mo)) for mo,co in pol.terms())
        adapted_norms.append(s.expand(poly))
    print('xy norm polynomials in adapted variables:')
    for n in adapted_norms:print(n)
    minor=s.Matrix([[s.Poly(adapted_norms[j],*adapted_vars).coeff_monomial(m) for j in [4,5]] for m in [F*V**2,G*V**2]])
    assert minor==s.eye(2)
    assert len(rref2(s.Matrix.hstack(*(mod(xy*q) for q in gens)))[0])==6
    orbits=s.Matrix.hstack(*[K*q for q in gens for K in [s.eye(24),Ii,Jj,Ii*Jj]])
    assert orbits.det()==3**6
    print('Last-two norm minor on (fv²,gv²):',minor.tolist())
    print('Six quaternionic-orbit determinant:',orbits.det())

    # Check the complete selected coefficient matrix printed in Appendix B.
    printed=s.Matrix([[1,0,0,1,0,0],[1,1,1,1,1,0],[1,1,1,0,0,1],[0,1,0,1,1,1],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    monomials=[F**3,F**2*G,F*G**2,G**3,F*V**2,G*V**2]
    certificate=s.Matrix([[s.Poly(poly,*adapted_vars).coeff_monomial(m) for poly in adapted_norms] for m in monomials])
    assert certificate==printed
    assert int(certificate.det())%2==1
    assert ker.cols==6 and ker[20:,:]==s.zeros(4,6)
    print('Complete printed six-by-six coefficient certificate verified.')

    # Exact integral identities displayed before reduction in Appendix B.
    common=b**3+c**3-d**3
    poly2=2*a**3+a**2*(2*b+2*c+d)+a*(3*b**2+3*c**2+4*c*d)+common+c**2*d+3*c*d**2
    poly3=a**2*(-b+2*c-2*d)+a*(-4*b*d+3*c**2-3*d**2)+common-3*b**2*d-b*d**2
    expected2=homotopy_cubic_column(poly2)
    expected2[20]=expected2[21]=expected2[22]=2
    expected3=homotopy_cubic_column(poly3)
    expected3[22]=2; expected3[23]=-2
    omega=(Ii-s.eye(24))*(Jj-s.eye(24))
    assert omega[:,21]==expected2
    assert omega[:,22]==expected3
    print('Both integral identities verified before reduction modulo two.')

    # Export only after the complete verification succeeds.
    if args.matrices is not None:
        with args.matrices.open("w", encoding="utf-8") as out:
            for name, matrix in [
                ("Hurewicz basis", B), ("i action", Ii),
                ("j action", Jj), ("mod2 fixed basis", ker),
            ]:
                out.write(name + "\n" + str(matrix.tolist()) + "\n\n")

    print('All checks passed: the degree-six lattice is L_A^6; boundary is injective.')


if __name__ == "__main__":
    main()
