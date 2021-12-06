from plane import Plane

class StressStateValues:

    def __init__(self, p=None, tau=None, mu_s=None, phi=None, sigma1=None, sigma2=None, sigma3=None):
        if sigma1 is not None and sigma2 is not None and sigma3 is not None:
            self.parse_from_sigmas(sigma1, sigma2, sigma3)
        elif p is not None and tau is not None and (mu_s is not None or phi is not None):
            self.parse_not_from_sigmas(p, tau, mu_s, phi)
        else:
            raise ValueError("You must set values either for (sigma1, sigma2, sigma3) or for (p, tau mus_s or phi)")

    def parse_from_sigmas(self, sigma1: float, sigma2: float, sigma3: float):
        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3
        self.p = -(sigma1 + sigma2 + sigma3) / 3
        self.tau = (sigma1 - sigma3) / 2
        self.mu_s = 2 * ((sigma2 - sigma3) / (sigma1 - sigma3)) - 1
        self.phi = (1 - self.mu_s) / 2

    def parse_not_from_sigmas(self, p, tau, mu_s, phi):
        mu_s, phi = self.parse_stress_shape(mu_s, phi)
        self.p = p
        self.tau = tau
        self.mu_s = mu_s
        self.phi = phi
        self.sigma1 = -p + tau * (1 - mu_s/3)
        self.sigma2 = -p + tau * (2 * mu_s/3)
        self.sigma3 = -p - tau * (1 + mu_s/3)

    def parse_stress_shape(self, mu_s, phi):
        if mu_s is None and phi is None:
            raise ValueError("Either phi or mu_s should be set to some value")
        if mu_s is None:
            mu_s = 1 - 2 * phi
            return mu_s, phi
        if phi is None:
            phi = (1 - mu_s)/2
            return mu_s, phi

    def __repr__(self):
        return (
            f"Values:\n"
            f"σ₁: {self.sigma1:.2f}\n"
            f"σ₂: {self.sigma2:.2f}\n"
            f"σ₃: {self.sigma3:.2f}\n"
            f"p: {self.p:.2f}\n"
            f"𝜏: {self.tau:.2f}\n"
            f"μ_σ: {self.mu_s:.2f}\n"
            f"φ: {self.phi:.2f}\n"
                )

class StressStateOrientation:

    def __init__(self, sigma1=None, sigma2=None, sigma3=None):

        def parse_single_sigma(maybe_sigma):
            if maybe_sigma is None:
                return None
            if isinstance(Plane, type(maybe_sigma)):
                return maybe_sigma
            if len(maybe_sigma) == 2:
                dr, dp = maybe_sigma
                return Plane(dr, dp)
            else:
                raise ValueError("Value must be either of type Plane or tuple of two float values")

        def return_third_axis_from_two(axis1, axis2):
            return axis1.get_perpendicular_between(axis2)

        sigma1 = parse_single_sigma(sigma1)
        sigma2 = parse_single_sigma(sigma2)
        sigma3 = parse_single_sigma(sigma3)

        sigmas = sigma1, sigma2, sigma3

        none_sigmas = sum(1 for i in sigmas if i is None)

        if none_sigmas > 1:
            raise ValueError ("Set orientations for at least two sigmas")

        if none_sigmas == 1:
            if sigma1 is None:
                sigma1 = return_third_axis_from_two(sigma2, sigma3)
            elif sigma2 is None:
                sigma2 = return_third_axis_from_two(sigma1, sigma3)
            elif sigma3 is None:
                sigma3 = return_third_axis_from_two(sigma1, sigma2)

        self.sigma1 = sigma1
        self.sigma2 = sigma2
        self.sigma3 = sigma3

    def __repr__(self):
        return (
            "Orientations:\n"
            f"σ₁: {self.sigma1.dir:.2f}° ∠{self.sigma1.dip:.2f}°\n"
            f"σ₂: {self.sigma2.dir:.2f}° ∠{self.sigma2.dip:.2f}°\n"
            f"σ₃: {self.sigma3.dir:.2f}° ∠{self.sigma3.dip:.2f}°"
                )


class StressState:
    values: StressStateValues
    orientation: StressStateOrientation

    def __init__(self, orientations, values):
        self.orientation = StressStateOrientation(**orientations)
        self.values = StressStateValues(**values)

    def __repr__(self):
        return (
            "STRESS STATE:\n\n"
            f"{self.orientation}\n"
            "\n"
            f"{self.values}"
        )
