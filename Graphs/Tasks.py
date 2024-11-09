from collections import deque

def bfs(graph, s, t, parent):
    visited = [False] * len(graph)
    queue = deque()
    queue.append(s)
    visited[s] = True

    while queue:
        u = queue.popleft()
        for v, capacity in enumerate(graph[u]):
            if not visited[v] and capacity > 0:
                queue.append(v)
                visited[v] = True
                parent[v] = u
                if v == t:
                    return True
    return False

def max_flow(graph, source, sink):
    parent = [-1] * len(graph)
    max_flow = 0

    while bfs(graph, source, sink, parent):
        path_flow = float("inf")
        s = sink
        while s != source:
            path_flow = min(path_flow, graph[parent[s]][s])
            s = parent[s]

        max_flow += path_flow
        v = sink
        while v != source:
            u = parent[v]
            graph[u][v] -= path_flow
            graph[v][u] += path_flow
            v = parent[v]

    return max_flow

def main():
    v, e = map(int, input().split())
    graph = [[0] * v for _ in range(v)]
    for _ in range(e):
        u, v, capacity = map(int, input().split())
        graph[u][v] += capacity
        graph[v][u] += capacity  # Учитываем обратные рёбра
    print(max_flow(graph, 0, v-1))

if __name__ == "__main__":
    main()
