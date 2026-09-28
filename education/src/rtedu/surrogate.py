"""Chapter 12: a small neural surrogate of R(T), in numpy, with hand-written
backprop, so that every line can be read. Input: (log T scaled, one-hot
group i); output: a softmax row over exit groups, trained by cross-entropy
against the tabulated rows. Whole physical states (temperatures) are held
out: no row of a held-out T is ever seen in training.
"""
import numpy as np


class MLP:
    def __init__(self, n_in, n_hidden, n_out, seed=0):
        rng = np.random.default_rng(seed)
        self.W1 = rng.normal(0, 1 / np.sqrt(n_in), (n_in, n_hidden)); self.b1 = np.zeros(n_hidden)
        self.W2 = rng.normal(0, 1 / np.sqrt(n_hidden), (n_hidden, n_out)); self.b2 = np.zeros(n_out)

    def forward(self, X):
        self.X = X; self.H = np.tanh(X @ self.W1 + self.b1)
        z = self.H @ self.W2 + self.b2
        z = z - z.max(axis=1, keepdims=True)
        self.P = np.exp(z); self.P /= self.P.sum(axis=1, keepdims=True)
        return self.P

    def loss(self, Y):
        return float(-(Y * np.log(self.P + 1e-12)).sum(axis=1).mean())

    def backward(self, Y, lr):
        n = Y.shape[0]
        dz = (self.P - Y) / n                      # d cross-entropy / d logits for softmax
        dW2 = self.H.T @ dz; db2 = dz.sum(axis=0)
        dH = dz @ self.W2.T * (1 - self.H ** 2)
        dW1 = self.X.T @ dH; db1 = dH.sum(axis=0)
        self.W1 -= lr * dW1; self.b1 -= lr * db1; self.W2 -= lr * dW2; self.b2 -= lr * db2

    def fit(self, X, Y, epochs=3000, lr=0.5):
        hist = []
        for _ in range(epochs):
            self.forward(X); hist.append(self.loss(Y)); self.backward(Y, lr)
        return np.array(hist)


def features(T, n_g, T_ref=4000.0):
    """One row of features per (T, group i): scaled log T and the one-hot group."""
    X = []
    for i in range(n_g):
        onehot = np.zeros(n_g); onehot[i] = 1.0
        X.append(np.concatenate([[np.log(T / T_ref)], onehot]))
    return np.array(X)


def dataset(Ts, Rs, n_g):
    X = np.vstack([features(T, n_g) for T in Ts]); Y = np.vstack([np.asarray(R) for R in Rs])
    return X, Y


def split_by_state(Ts, held_out):
    """Indices of training and held-out temperatures: entire states withheld."""
    held = [i for i, T in enumerate(Ts) if T in held_out]
    train = [i for i in range(len(Ts)) if i not in held]
    return train, held


def predict_R(model, T, n_g):
    return model.forward(features(T, n_g))
