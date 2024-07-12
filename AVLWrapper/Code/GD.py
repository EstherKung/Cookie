import numpy as np
import matplotlib.pyplot as plt

def ms_error(y_true, y_pred):
    cost = np.sum((y_true - y_pred)**2) / len(y_true)
    return cost

def gradient_descent(x, y, z, iterations = 1000, learning_rate = 0.0001, threshold = 1e-6):
    w = 0.1
    b = 0.01
    iter = iterations
    lr = learning_rate
    n = float(len(x))

    costs = []; weights = []
    cost_prev = None
    
    for i in range(iter):
        y_pred = w*x + b
        cost = ms_error(y, y_pred)

        if cost_prev and abs(cost_prev - cost) <= threshold:
            break

        cost_prev = cost

        costs.append(cost)
        weights.append(w)

        dJdw = -(2/n) * sum(x * (y-y_pred))
        dJdb = -(2/n) * sum(y-y_pred)

        w = w - (lr * dJdw)
        b = b - (lr * dJdb)

        print(f"Iteration{i+1}: Cost {cost}, Weight {w}, Bias {b}")

    return w, b

#x, input, elevator size
#y, output, elevator deflection
#z, output, static margin