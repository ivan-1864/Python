import numpy as np

class LinearRegression:
    def __init__(self, **kwargs):
        self.coef_ = None
        pass

    def fit(self, x: np.array, y: np.array):
        X = np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)
        self.coef_ = np.linalg.inv(X.T @ X) @ X.T @ y 
        pass

    def predict(self, x: np.array):
        X = np.concatenate((x, np.ones((x.shape[0], 1))), axis=1)
        return X @ self.coef_

def r2(y_true, y_pred):
  De = np.sum(np.power(y_true - y_pred, 2))
  Dz = np.sum(np.power(y_true - y_true.mean(), 2))
  return 1 - De / Dz
x=[]
reg=LinearRegression()
x.append(np.load('1.npy'))
x.append(np.load('2.npy'))
x.append(np.load('3.npy'))
x.append(np.load('4.npy'))
x.append(np.load('5.npy'))
for i in range(5):
    X = x[i][:, 0].reshape((-1, 1))
    y=x[i][:, 1]
    reg.fit(X,y)
    print(r2(y, reg.predict(X)))
