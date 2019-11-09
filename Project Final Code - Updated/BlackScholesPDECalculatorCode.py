from math import exp, floor
import matrixOperations


class Options(object):
    def __init__(self, S, K, r, q, T, sigma):
        self.S = float(S)
        self.K = float(K)
        self.r = float(r)
        self.q = float(q)
        self.T = float(T)
        self.sigma = float(sigma)


class EuropeanOptions(Options):
    def ExplicitFDM(self, M, N):
        # 1 Compute deltat = T/n; deltaS = Smax/m where Smax = 2*K
        t_delta = self.T / N
        S_max = 2 * self.K
        S_delta = S_max / M
        # 2 Compute (Call option)fn,j=max(j*deltaS-K,0);
        #  (Put option)fn,j=max(K-j*deltaS),for j=0,1,...,M
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - self.K, 0)
            Fp[j][N] = max(self.K - j * S_delta, 0)
        # 3 For i=N-1, N-2,...,1,0, repeat 3.1 and 3.2
        for i in range(N - 1, -1, -1):
            # 3.1 Compute vector Fi_hat=A*Fi+1,where:
            for j in range(1, M):
                a = 0.5 * t_delta * (self.sigma * self.sigma * j * j
                                     - (self.r - self.q) * j)
                b = 1 - t_delta * (self.sigma * self.sigma * j * j + self.r)
                c = 0.5 * t_delta * (self.sigma * self.sigma * j * j
                                     + (self.r - self.q) * j)
                Fc[j][i] = (a * Fc[j - 1][i + 1]
                            + b * Fc[j][i + 1]
                            + c * Fc[j + 1][i + 1])
                Fp[j][i] = (a * Fp[j - 1][i + 1]
                            + b * Fp[j][i + 1]
                            + c * Fp[j + 1][i + 1])
                # 3.2 Compute vector Fi
            Fc[0][i] = 0
            Fc[M][i] = S_max - self.K * exp(-self.r * (N - i) * t_delta)
            Fp[0][i] = self.K * exp(-self.r * (N - i) * t_delta)
            Fp[M][i] = 0
        # 4 Find k, such that k*(deltaS)<=S<=(k+1)*(deltaS), i.e. k=[S/deltaS]
        #  [.] represents the floor of .
        k = int(floor(self.S / S_delta))
        # 5 Option price: V=f0,k + ((f0,k+1 - f0,k)/deltaS)*(S-k*deltaS)
        Vc = (Fc[k][0] + (Fc[k + 1][0] - Fc[k][0])
              / S_delta * (self.S - k * S_delta))
        Vp = (Fp[k][0] + (Fp[k + 1][0] - Fp[k][0])
              / S_delta * (self.S - k * S_delta))
        return {'call': Vc, 'put': Vp}

    def ImplicitFDM(self, M, N):
        # 1
        t_delta = self.T / N
        S_max = 2 * self.K
        S_delta = S_max / M
        # 2
        Fc_hat = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp_hat = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - self.K, 0)
            Fp[j][N] = max(self.K - j * S_delta, 0)
        # 3
        A = [[0 for _ in range(M + 1)] for _ in range(M + 1)]
        A[0][0] = 1
        for j in range(1, M):
            A[j][j - 1] = 0.5 * t_delta * ((self.r - self.q) * j
                                           - self.sigma * self.sigma * j * j)
            A[j][j] = 1 + t_delta * (self.sigma * self.sigma
                                     * j * j + self.r)
            A[j][j + 1] = -0.5 * t_delta * (self.sigma * self.sigma
                                            * j * j + (self.r - self.q) * j)
        A[M][M] = 1
        A_inv = matrixOperations.inverse(A)

        for i in range(N - 1, -1, -1):
            # 3.1
            for j in range(1, M):
                Fc_hat[j][i + 1] = Fc[j][i + 1]
                Fp_hat[j][i + 1] = Fp[j][i + 1]
            Fc_hat[0][i + 1] = 0
            Fc_hat[M][i + 1] = S_max - self.K * exp(-self.r * (N - i) * t_delta)
            Fp_hat[0][i + 1] = self.K * exp(-self.r * (N - i) * t_delta)
            Fp_hat[M][i + 1] = 0
            # 3.2

            Fc_hat_col = [[Fc_hat[j][i + 1]] for j in range(M + 1)]
            Fp_hat_col = [[Fp_hat[j][i + 1]] for j in range(M + 1)]
            Fc_col = matrixOperations.multiply(A_inv, Fc_hat_col)
            Fp_col = matrixOperations.multiply(A_inv, Fp_hat_col)
            for j in range(M + 1):
                Fc[j][i] = Fc_col[j][0]
                Fp[j][i] = Fp_col[j][0]
        # 4
        k = int(floor(self.S / S_delta))
        # 5
        Vc = (Fc[k][0] + (Fc[k + 1][0] - Fc[k][0])
              / S_delta * (self.S - k * S_delta))
        Vp = (Fp[k][0] + (Fp[k + 1][0] - Fp[k][0])
              / S_delta * (self.S - k * S_delta))
        return {"call": Vc, "put": Vp}

    def CrankNicolsonM(self, M, N):
        # 1
        t_delta = self.T / N
        S_max = 2 * self.K
        S_delta = S_max / M
        # 2
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - self.K, 0)
            Fp[j][N] = max(self.K - j * S_delta, 0)

        # 3
        M1 = [[0 for _ in range(M + 1)] for _ in range(M + 1)]
        M1[0][0] = 1
        M1[M][M] = 1

        for j in range(1, M):
            M1[j][j - 1] = -0.25 * t_delta * (self.sigma * self.sigma
                                              * j * j - (self.r - self.q) * j)
            M1[j][j] = 1 + 0.5 * t_delta * (self.sigma * self.sigma * j * j + self.r)
            M1[j][j + 1] = -0.25 * t_delta * (self.sigma * self.sigma
                                              * j * j + (self.r - self.q) * j)
        M1_inv = matrixOperations.inverse(M1)

        for i in range(N - 1, -1, -1):
            # 3.1 and 3.2:
            bc = [[0] for _ in range(M + 1)]
            bp = [[0] for _ in range(M + 1)]
            for j in range(1, M):
                alpha = 0.25 * t_delta * (self.sigma * self.sigma * j * j
                                          - (self.r - self.q) * j)
                beta = -0.5 * t_delta * (self.sigma * self.sigma * j * j + self.r)
                gamma = 0.25 * t_delta * (self.sigma * self.sigma * j * j
                                          + (self.r - self.q) * j)
                bc[j][0] = (alpha * Fc[j - 1][i + 1]
                            + (1 + beta) * Fc[j][i + 1]
                            + gamma * Fc[j + 1][i + 1])
                bp[j][0] = (alpha * Fp[j - 1][i + 1]
                            + (1 + beta) * Fp[j][i + 1]
                            + gamma * Fp[j + 1][i + 1])
            bc[0][0] = 0
            bc[M][0] = S_max - self.K * exp(-self.r * (N - i) * t_delta)
            bp[0][0] = self.K * exp(-self.r * (N - i) * t_delta)
            bp[M][0] = 0

            # 3.3
            Fc_col = matrixOperations.multiply(M1_inv, bc)
            Fp_col = matrixOperations.multiply(M1_inv, bp)
            for j in range(M + 1):
                Fc[j][i] = Fc_col[j][0]
                Fp[j][i] = Fp_col[j][0]
        # 4
        k = int(floor(self.S / S_delta))
        # 5
        Vc = (Fc[k][0] + (Fc[k + 1][0] - Fc[k][0])
              / S_delta * (self.S - k * S_delta))
        Vp = (Fp[k][0] + (Fp[k + 1][0] - Fp[k][0])
              / S_delta * (self.S - k * S_delta))
        return {"call": Vc, "put": Vp}
