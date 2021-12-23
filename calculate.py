from math import tan, sin, cos, radians, atan2, degrees
from dataclasses import dataclass

from rich.progress import track

from plane import Plane



@dataclass
class StressOnPlane:
    plane: Plane
    s_nn_reduced: float
    s_nm_reduced: float
    s_nt_reduced: float
    tau: float | None = None
    p: float | None = None

    @property
    def s_nn(self):
        if self.tau is None or self.p is None:
            raise ValueError("Values for 'p' and 'tau' should be set to calculate real stress")
        return -self.p + self.tau * self.s_nn_reduced

    @property
    def s_nm(self):
        if self.tau is None or self.p is None:
            raise ValueError("Values for 'p' and 'tau' should be set to calculate real stress")
        return self.tau * self.s_nm_reduced

    @property
    def s_nt(self):
        if self.tau is None or self.p is None:
            raise ValueError("Values for 'p' and 'tau' should be set to calculate real stress")
        return self.tau * self.s_nt_reduced

    @property
    def phi(self):
        return degrees(atan2(self.s_nt_reduced, self.s_nm_reduced))# + 90

    @property
    def tau_n(self):
        return (self.s_nt**2 + self.s_nm**2)**.5

    @property
    def tau_n_reduced(self):
        return (self.s_nt_reduced**2 + self.s_nm_reduced**2)**.5


def rotate_plane(s1, s3, p):
    """Return new plane rotated so s1 has dr=0, dp=0, s3 has dr=90, dp=0."""
    p_new = Plane(p.dir, p.dip)
    s1_new = Plane(s1.dir, s1.dip)

    s2 = s1.get_perpendicular_between(s3)
    s2_n = s2.normal()

    dr = s2_n.dir
    dp = s2_n.dip

    p_new.rotate_cos2(-dr)
    p_new.rotate_cos3(-dp)

    s1_new.rotate_cos2(-dr)
    s1_new.rotate_cos3(-dp)

    p_new.rotate_cos2(-s1_new.dir)

    return p_new

def end_of_direction(plane, phi, lenght=5):
    alpha = plane.dir
    beta = plane.dip

    end_plane = Plane(phi, 90-lenght)

    end_plane.rotate_cos3(beta-90)
    end_plane.rotate_cos2(alpha)

    return end_plane

def chronick(i, n):
    return 1 if n==i else 0

def s_ni_reduced(i, n, mu_s, l_i1, l_i3, l_n1, l_n3):
    return (1-mu_s) * l_n1 * l_i1 - (1+mu_s) * l_n3 * l_i3 + chronick(i, n) * mu_s * (2/3.0)

def get_s_tau(plane, mu_s, phi=0):
    plane_m = Plane(plane.dir, plane.dip+90).rotated(phi)
    plane_t = Plane(plane.dir, plane.dip+90).rotated(phi + 90)

    l_n1 = plane.cos1
    l_n3 = plane.cos3

    l_m1 = plane_m.cos1
    l_m3 = plane_m.cos3

    l_t1 = plane_t.cos1
    l_t3 = plane_t.cos3

    s_nn_reduced = s_ni_reduced(i=0,n=0, mu_s=mu_s, l_i1=l_n1, l_i3=l_n3, l_n1=l_n1, l_n3=l_n3)
    s_nm_reduced = s_ni_reduced(i=0,n=1, mu_s=mu_s, l_i1=l_m1, l_i3=l_m3, l_n1=l_n1, l_n3=l_n3)
    s_nt_reduced = s_ni_reduced(i=0,n=1, mu_s=mu_s, l_i1=l_t1, l_i3=l_t3, l_n1=l_n1, l_n3=l_n3)

    return s_nn_reduced, s_nm_reduced, s_nt_reduced

def calculate_stress(plane, stress_state):
    plane_rotated = rotate_plane(
            stress_state.orientation.sigma1,
            stress_state.orientation.sigma3,
            plane,
        )

    s_nn_reduced, s_nm_reduced, s_nt_reduced = get_s_tau(
            plane=plane_rotated,
            mu_s=stress_state.values.mu_s,
            )

    return StressOnPlane(
        plane=plane,
        s_nn_reduced=s_nn_reduced,
        s_nm_reduced=s_nm_reduced,
        s_nt_reduced=s_nt_reduced,
        tau=stress_state.values.tau,
        p=stress_state.values.p,
        )


def calculate_stress_on_planes(stress_state):
    stresses_on_plane = []

    for dr in track(range(0, 361, 2)):
        for dp in range(0, 91, 1):
            stresses_on_plane.append(calculate_stress(Plane(dr, dp), stress_state))

    return stresses_on_plane

def calculate_stresses_on_fractures(stress_state, fractures):
    stresses_on_fractures = []

    for fracture in fractures:
        stresses_on_fractures.append(calculate_stress(fracture.normal(), stress_state))

    return stresses_on_fractures

def fracture_criteria_reduced(stress_on_plane, tau_f, k_f):
    s_nn = stress_on_plane.s_nn_reduced
    tau_n = stress_on_plane.tau_n_reduced
    tau2 = tau_f  - k_f * s_nn

    return 0 if tau_n > tau2 else 1
