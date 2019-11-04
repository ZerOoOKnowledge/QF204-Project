import math
import matrixOperations


class Solver:
    def __init__(self, S, K, r, q, T, sigma):
        self.S = S
        self.K = K
        self.r = r
        self.q = q
        self.T = T
        self.sigma = sigma

    def explicit(self, N, M):
        S = self.S
        K = self.K
        r = self.r
        q = self.q
        T = self.T
        sigma = self.sigma

        # Step 1:
        t_delta = T / N
        S_max = 2 * K
        S_delta = S_max / M

        # Step 2:
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - K, 0)
            Fp[j][N] = max(K - j * S_delta, 0)

        # Step 3:
        for i in range(N - 1, -1, -1):
            # Step 3.1:
            for j in range(1, M):
                a = 0.5 * t_delta * (sigma * sigma * j * j - (r - q) * j)
                b = 1 - t_delta * (sigma * sigma * j * j + r)
                c = 0.5 * t_delta * (sigma * sigma * j * j + (r - q) * j)
                Fc[j][i] = a * Fc[j - 1][i + 1] + b * Fc[j][i + 1] + c * Fc[j + 1][i + 1]
                Fp[j][i] = a * Fp[j - 1][i + 1] + b * Fp[j][i + 1] + c * Fp[j + 1][i + 1]

            # Step 3.2:
            Fc[0][i] = 0
            Fc[M][i] = S_max - K * math.exp(-r * (N - i) * t_delta)
            Fp[0][i] = K * math.exp(-r * (N - i) * t_delta)
            Fp[M][i] = 0

        # Step 4:
        k = int(math.floor(S / S_delta))

        # Step 5:
        Vc = Fc[k][0] + (Fc[k + 1][0] - Fc[k][0]) / S_delta * (S - k * S_delta)
        Vp = Fp[k][0] + (Fp[k + 1][0] - Fp[k][0]) / S_delta * (S - k * S_delta)

        return {"call": Vc, "put": Vp}

    def implicit(self, N, M):
        S = self.S
        K = self.K
        r = self.r
        q = self.q
        T = self.T
        sigma = self.sigma

        # Step 1:
        t_delta = T / N
        S_max = 2 * K
        S_delta = S_max / M

        # Step 2:
        Fc_hat = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp_hat = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - K, 0)
            Fp[j][N] = max(K - j * S_delta, 0)

        # Step 3:
        for i in range(N - 1, -1, -1):
            # Step 3.1:
            for j in range(1, M):
                Fc_hat[j][i + 1] = Fc[j][i + 1]
                Fp_hat[j][i + 1] = Fp[j][i + 1]
            Fc_hat[0][i + 1] = 0
            Fc_hat[M][i + 1] = S_max - K * math.exp(-r * (N - i) * t_delta)
            Fp_hat[0][i + 1] = K * math.exp(-r * (N - i) * t_delta)
            Fp_hat[M][i + 1] = 0

            # Step 3.2:
            A = [[0 for _ in range(M + 1)] for _ in range(M + 1)]
            A[0][0] = 1
            for j in range(1, M):
                A[j][j - 1] = 0.5 * t_delta * ((r - q) * j - sigma * sigma * j * j)
                A[j][j] = 1 + t_delta * (sigma * sigma * j * j + r)
                A[j][j + 1] = -0.5 * t_delta * (sigma * sigma * j * j + (r - q) * j)
            A[M][M] = 1
            A_inv = matrixOperations.inverse(A)
            Fc_hat_col = [[Fc_hat[j][i + 1]] for j in range(M + 1)]
            Fp_hat_col = [[Fp_hat[j][i + 1]] for j in range(M + 1)]
            Fc_col = matrixOperations.multiply(A_inv, Fc_hat_col)
            Fp_col = matrixOperations.multiply(A_inv, Fp_hat_col)
            for j in range(M + 1):
                Fc[j][i] = Fc_col[j][0]
                Fp[j][i] = Fp_col[j][0]

        # Step 4:
        k = int(math.floor(S / S_delta))

        # Step 5:
        Vc = Fc[k][0] + (Fc[k + 1][0] - Fc[k][0]) / S_delta * (S - k * S_delta)
        Vp = Fp[k][0] + (Fp[k + 1][0] - Fp[k][0]) / S_delta * (S - k * S_delta)

        return {"call": Vc, "put": Vp}

    def crankNicolson(self, N, M):
        S = self.S
        K = self.K
        r = self.r
        q = self.q
        T = self.T
        sigma = self.sigma

        # Step 1:
        t_delta = T / N
        S_max = 2 * K
        S_delta = S_max / M

        # Step 2:
        Fc = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        Fp = [[0 for _ in range(N + 1)] for _ in range(M + 1)]
        for j in range(M + 1):
            Fc[j][N] = max(j * S_delta - K, 0)
            Fp[j][N] = max(K - j * S_delta, 0)

        # Step 3:
        for i in range(N - 1, -1, -1):
            # Step 3.1, 3.2:
            bc = [[0] for _ in range(M + 1)]
            bp = [[0] for _ in range(M + 1)]
            for j in range(1, M):
                alpha = 0.25 * t_delta * (sigma * sigma * j * j - (r - q) * j)
                beta = -0.5 * t_delta * (sigma * sigma + r)
                gamma = 0.25 * t_delta * (sigma * sigma * j * j + (r - q) * j)
                bc[j][0] = alpha * Fc[j - 1][i + 1] + beta * Fc[j][i + 1] + gamma * Fc[j + 1][i + 1]
                bp[j][0] = alpha * Fp[j - 1][i + 1] + beta * Fp[j][i + 1] + gamma * Fp[j + 1][i + 1]
            bc[0][0] = 0
            bc[M][0] = S_max - K * math.exp(-r * (N - i) * t_delta)
            bp[0][0] = K * math.exp(-r * (N - i) * t_delta)
            bp[M][0] = 0

            # Step 3.3:
            M1 = [[0 for _ in range(M + 1)] for _ in range(M + 1)]
            M1[0][0] = 1
            M1[M][M] = 1

            for j in range(1, M):
                M1[j][j - 1] = 0.25 * t_delta * (sigma * sigma * j * j - (r - q) * j)
                M1[j][j] = -0.5 * t_delta * (sigma * sigma + r)
                M1[j][j + 1] = 0.25 * t_delta * (sigma * sigma * j * j + (r - q) * j)

            M1_inv = matrixOperations.inverse(M1)
            Fc_col = matrixOperations.multiply(M1_inv, bc)
            Fp_col = matrixOperations.multiply(M1_inv, bp)
            for j in range(M + 1):
                Fc[j][i] = Fc_col[j][0]
                Fp[j][i] = Fp_col[j][0]

        # Step 4:
        k = int(math.floor(S / S_delta))

        # Step 5:
        Vc = Fc[k][0] + (Fc[k + 1][0] - Fc[k][0]) / S_delta * (S - k * S_delta)
        Vp = Fp[k][0] + (Fp[k + 1][0] - Fp[k][0]) / S_delta * (S - k * S_delta)

        return {"call": Vc, "put": Vp}
