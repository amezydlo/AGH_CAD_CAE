import numpy as np
import random

np.set_printoptions(threshold=np.inf)


def generate_empty_maze(n):
    """Generates an empty n x n matrix filled with walls (1)."""
    maze = np.ones((n, n), dtype=int)
    return maze


def find_neighbors(x, y, maze):
    """Returns a list of neighbors that are two cells away from the current position."""
    neighbors = []
    if x > 1 and maze[x - 2, y] == 1:
        neighbors.append((x - 2, y))  # Up
    if x < maze.shape[0] - 2 and maze[x + 2, y] == 1:
        neighbors.append((x + 2, y))  # Down
    if y > 1 and maze[x, y - 2] == 1:
        neighbors.append((x, y - 2))  # Left
    if y < maze.shape[1] - 2 and maze[x, y + 2] == 1:
        neighbors.append((x, y + 2))  # Right
    return neighbors


def generate_maze_dfs(n):
    """Generates a random n x n maze using the DFS algorithm."""
    maze = generate_empty_maze(n)

    start_x, start_y = (random.randrange(1, n, 2), random.randrange(1, n, 2))
    maze[start_x, start_y] = 0  # Mark the starting point as a path (0)

    stack = [(start_x, start_y)]

    while stack:
        x, y = stack[-1]
        neighbors = find_neighbors(x, y, maze)

        if neighbors:
            new_x, new_y = random.choice(neighbors)
            maze[(x + new_x) // 2, (y + new_y) // 2] = 0
            maze[new_x, new_y] = 0
            stack.append((new_x, new_y))
        else:
            stack.pop()

    return maze


def generate_knot_vectors(n, k):
    """Generates a knot vector for a given axis."""
    knot_vector = np.linspace(0, n - 2, k, dtype=int)
    return knot_vector


# 2 * powtorzenia + reszta = 36 + 1 + 2
# 39 - 6 = 33

# liczba spline = 10 = knot_vec - 1 - 2
# 13 = knot_vec
# 13 - 6 = 7

# Parameters
n = 46  # Maze size (must be odd)

# 69 + 6 - 1 - 2 = 72 / 2 = 36

# 36 = len_k - 2 - 1 -> 39

# 40 - 2 - 1 = 37
# Maze generation
maze = generate_maze_dfs(n)

print("Maze:")
print(maze)
print("Knot vectors")

# print(np.linspace(0, 20, ))
print([0] * 3 + np.linspace(1, n-3 , n-3, dtype=int).tolist() + [n - 2] * 3)

