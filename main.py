import networkx as nx
import matplotlib.pyplot as plt
import random
import math
EPS=10e-5
class KamadaKawai:
    def __init__(self, graph, K=30, L0=30, INF=1000):
        self.n = graph.number_of_nodes()
        self.graph = graph
        self.INF = INF
        self.x = [random.randint(0, 200) for i in range(self.n)]
        self.y = [random.randint(0, 200) for j in range(self.n)]
        self.delta = [[0, 0] for i in range(self.n)]
        # 最短経路問題
        self.warshall_floyd()
    
        self.l = [[ 0 if self.dist[i][j] == 0 else L0 *  self.dist[i][j]       for i in range(self.n)] for j in range(self.n)]
        self.k = [[ 0 if self.dist[i][j] == 0 else K /  (self.dist[i][j] ** 2) for i in range(self.n)] for j in range(self.n)]

    def warshall_floyd(self):
        self.dist = [[self.INF for i in range(self.n)] for j in range(self.n)]
        for s, t in graph.edges():
            self.dist[s][t] = self.dist[t][s] = 1
        for i in range(self.n):
            self.dist[i][i] = 0
        for k in range(self.n):
            for i in range(self.n):
                for j in range(self.n):
                    self.dist[i][j] = min(self.dist[i][j], self.dist[i][k] + self.dist[k][j])

    def norm(self, i, j):
        square = (self.x[i] - self.x[j]) ** 2 + \
                 (self.y[i] - self.y[j]) ** 2
        return math.sqrt(square)

    def determine_point(self):
        # Eをすべての点で微分し、一番大きくなる点を選ぶ
        max_point = 0
        max_delta = 0
        for i in range(self.n):
            ex_pi = 0
            ey_pi = 0
            for j in range(self.n):
                if (i == j):
                    continue
                norm = self.norm(i, j)
                # Eのxの微分
                ex_pi += self.k[i][j] * (self.x[i] - self.x[j]) * (1 - self.l[i][j] / norm)
                # Eのyの微分
                ey_pi += self.k[i][j] * (self.y[i] - self.y[j]) * (1 - self.l[i][j] / norm)
            self.delta[i] = [ex_pi, ey_pi]
            delta = math.sqrt(ex_pi ** 2 + ey_pi ** 2)
            if max_delta < delta:
                max_delta = delta
                max_point = i
        return max_point
    
    def calc_differential(self, op_p):
        exx = exy = eyy = 0 
        for i in range(self.n):
            if (i == op_p):
                continue
            norm = self.norm(i, op_p)
            dif_x = self.x[op_p] - self.x[i]
            dif_y = self.y[op_p] - self.y[i]
            ex, ey = self.delta[op_p]
            exy += self.k[op_p][i] *        self.l[op_p][i] * dif_x * dif_y / pow(norm, 3)
            exx += self.k[op_p][i] * (1.0 - self.l[op_p][i] * dif_y * dif_y / pow(norm, 3))
            eyy += self.k[op_p][i] * (1.0 - self.l[op_p][i] * dif_x * dif_x / pow(norm, 3))
        return ex, ey, exx, exy, eyy


graph = nx.karate_club_graph()
kk = KamadaKawai(graph)
# 最適化する点
op_i = kk.determine_point() 
# 一個前のループの差分
ex_pi, ey_pi = kk.delta[op_i]
delta = math.sqrt(ex_pi ** 2 + ey_pi ** 2)
while delta > EPS:
    ex, ey, exx, exy, eyy = kk.calc_differential(op_i)
    d = exx * eyy - exy * exy
    dx = - ( eyy * ex - exy * ey) / d
    dy = - (-exy * ex + exx * ey) / d
    kk.x[op_i] += dx
    kk.y[op_i] += dy
    delta = math.sqrt(ex**2 + ey**2)
    print(delta)
    op_i = kk.determine_point()
pos = {i: (kk.x[i], kk.y[i]) for i in range(kk.n)}
nx.draw(kk.graph, pos=pos, with_labels = True)
plt.show()