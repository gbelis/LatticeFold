# Contains functions to manipulate structure object on a lattice
import numpy as np



rotation_matrix = lambda theta: (np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]]).astype(int))

def rotation(shape,center=np.array([0,0]).reshape((2,1)),theta=0):
    """
    Apply a 2D rotation to a shape around a specified center.

    Parameters:
        shape (np.ndarray): A 2D array representing the coordinates of the shape.
        center (np.ndarray): A 2x1 array specifying the center of rotation.
        theta (float): The angle of rotation in radians.

    Returns:
        np.ndarray: The rotated shape as a 2D array.
    """
    return np.dot(rotation_matrix(theta),shape-center.reshape((2,1)))+center.reshape((2,1))

def are_neighbors(p1, p2) -> bool:
    """
    Check if two points are neighbors on a 2D lattice.

    Parameters:
        p1 (tuple): Coordinates of the first point.
        p2 (tuple): Coordinates of the second point.

    Returns:
        bool: True if the points are neighbors, False otherwise.
    """
    x1, y1 = p1
    x2, y2 = p2
    return abs(x1 - x2) + abs(y1 - y2) == 1

def get_neighbors(position):
    """
    Get the coordinates of all neighbors of a given position on a 2D lattice.

    Parameters:
        position (tuple): The (x, y) coordinates of the point.

    Returns:
        list of tuple: A list of tuples representing the neighboring positions.
    """
    x, y = position
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

def sterical_clash(structure) -> bool:
    """
    Check if there is a sterical clash in a given set of positions.

    Parameters:
        pos (np.ndarray): A 2D array of positions.

    Returns:
        bool: True if a sterical clash is detected, False otherwise.
    """
    return not len(np.unique(structure, axis=0))==len(structure)

def chain_connected(strucutre):
    """
    Check if a chain of points in a 2D lattice is connected.

    Parameters:
        structure (np.ndarray or list): Array/list of points representing the chain.

    Returns:
        bool: True if the chain is connected, False otherwise.
    """
    # Iterate through the chain and verify connectivity between consecutive points
    for i in range(len(strucutre) - 1):
        if not are_neighbors(strucutre[i], strucutre[i + 1]):
            return False
    return True


def generate_connected_points(n, start=(0, 0)):   ### NEEDS IMPROVING. NOT SURE CONCEPTUALY HOW TO SOLVE THIS.
    """
    Generates a random sequence of connected points in 2D without repeating points,
    using backtracking to avoid getting stuck.
    
    :param n: Number of points in the sequence
    :return: List of (x, y) coordinates
    """
    sequence = [start]
    
    while len(sequence) < n:
        last_x, last_y = sequence[-1]
        moves = [(1, 0), (-1, 0), (0, 1), (0, -1)]
        
        for dx, dy in np.random.permutation(moves):
            new_point = (last_x + dx, last_y + dy)
            if new_point not in sequence:
                sequence.append(new_point)
                break
        else:
            # No valid moves found, backtrack
            sequence.pop()
    
    return np.array(sequence)