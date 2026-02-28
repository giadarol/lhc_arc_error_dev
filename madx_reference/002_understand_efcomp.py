from math import factorial

from cpymad.madx import Madx


mad = Madx()

k0l = 0.01
r_ref = 0.3
dk2r = 0.2
dk2sr = 0.3
dk3r = 0.4
dk3sr = 0.5

src = (f"""
k0l={k0l};
r_ref={r_ref};
dk2r={dk2r};
dk2sr={dk2sr};
dk3r={dk3r};
dk3sr={dk3sr};
"""
+
"""

elm: multipole, knl:={k0l}, ksl={0.};

seq: sequence, l=1;
elm, at=0;
endsequence;

beam;
use,sequence=seq;

select,pattern=elm,flag=error;
efcomp,order=0, radius=r_ref,
    ! dkn={0.01,-0.02,0.03,0},
    ! dks={-0.01,0.02,-0.03,0.3,5},
    dknr={0.0,-0.04,dk2r,dk3r},
    dksr={-0.03,0.05,dk2sr,dk3sr};
""")

mad.input(src)

dkn_mad = mad.sequence.seq.expanded_elements['elm'].field_errors.dkn
dks_mad = mad.sequence.seq.expanded_elements['elm'].field_errors.dks

expected_k2l = dk2r * k0l * r_ref**(0 - 2) * factorial(2) / factorial(0)
print(f"Expected k2l: {expected_k2l:.7e}")
print(f"   Found k2l: {dkn_mad[2]:.7e}")

expected_k2sl = dk2sr * k0l * r_ref**(0 - 2) * factorial(2) / factorial(0)
print(f"Expected k2sl: {expected_k2sl:.7e}")
print(f"   Found k2sl: {dks_mad[2]:.7e}")

expected_k3l = dk3r * k0l * r_ref**(0 - 3) * factorial(3) / factorial(0)
print(f"Expected k3l: {expected_k3l:.7e}")
print(f"   Found k3l: {dkn_mad[3]:.7e}")

expected_k3sl = dk3sr * k0l * r_ref**(0 - 3) * factorial(3) / factorial(0)
print(f"Expected k3sl: {expected_k3sl:.7e}")
print(f"   Found k3sl: {dks_mad[3]:.7e}")