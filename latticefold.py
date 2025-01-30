import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

# Define constants
ENERGY_HH = -1  # Energy for H-H contact
GRID_SIZE = 20  # Size of the lattice grid


# Utils
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

def sterical_clash(pos) -> bool:
    """
    Check if there is a sterical clash in a given set of positions.

    Parameters:
        pos (np.ndarray): A 2D array of positions.

    Returns:
        bool: True if a sterical clash is detected, False otherwise.
    """
    return not len(np.unique(pos, axis=0))==len(pos)

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




class Protein:
    def __init__(self,sequence,initial_structure=None) -> None:
        self.sequence=sequence
        if initial_structure is None:
            self.structure=generate_connected_points(len(sequence))
        else: self.structure=initial_structure

    
    def get_moves(self, position, check=True, rules='all'):
        """
        Generate possible moves for a specific position in the current protein structure.

        Parameters:
            position (int): The index of the residue to generate moves for.
            check (bool): Whether to validate moves for chain connectivity and sterical clash.
            rules (str or list(str)): The rules for generating moves (combination of: 'small', 'rotation', 'crankshaft', or 'all').

        Returns:
            list of np.ndarray: A list of possible moves as 2D numpy arrays.
        """
        moves = []

        if 'rotation' in rules or rules == 'all':
            if position > 0:
                pos = self.structure[position]
                for theta in np.arange(1, 4) / 2 * np.pi:
                    move = self.structure.T.copy()
                    move[:, position:] = rotation(move[:, position:], pos.reshape(2, 1), theta)
                    if not check:
                        moves.append(move.T)
                    elif chain_connected(move.T) and not sterical_clash(move.T):
                        moves.append(move.T)

        if 'small' in rules or rules == 'all':
            if position == 0 or position == len(self.structure) - 1:  # End points
                neighbors = get_neighbors(self.structure[position])
                for nb in neighbors:
                    if not (nb == self.structure[position]).all():
                        move = self.structure.copy()
                        move[position] = nb
                        if not check:
                            moves.append(move)
                        elif chain_connected(move) and not sterical_clash(move):
                            moves.append(move.copy())
            elif 1 <= position < len(self.structure) - 1:  # Internal positions
                move = self.structure.copy()
                if np.abs(move[position - 1] - move[position + 1])[0] == 1:
                    move[position] = move[position - 1] - move[position] + move[position + 1]
                    if not check:
                        moves.append(move)
                    elif chain_connected(move) and not sterical_clash(move):
                        moves.append(move)

        if 'crankshaft' in rules or rules == 'all':
            if 1 <= position < len(self.structure) - 2:  # Ensure valid crankshaft rotation range
                move = self.structure.copy()
                v1, v2 = move[position - 1], move[position + 2]
                if np.linalg.norm(v1 - v2) == 2:  # Ensuring valid crankshaft rotation
                    move[position], move[position + 1] = move[position + 1], move[position]  # Swap positions
                    if not check:
                        moves.append(move)
                    elif chain_connected(move) and not sterical_clash(move):
                        moves.append(move.copy())

        return moves


    def get_energy(self, structure=None):
        """
        Calculate the total energy of the protein structure based on H-H contacts.

        Parameters:
            structure (np.ndarray, optional): The structure for which to calculate energy. Defaults to the current structure.

        Returns:
            int: The total energy of the structure.
        """
        energy = 0
        if structure is None:
            structure = self.structure  # Default to the current structure
        
        for i, pos in enumerate(structure):
            if self.sequence[i] == 'H':  # Only calculate for H residues
                neighbors = get_neighbors(pos)
                for neighbor in neighbors:
                    if any((neighbor == structure).all(axis=1)):  # Check if neighbor exists in structure
                        j = np.where((structure == neighbor).all(axis=1))[0][0]  # Get index of the neighbor
                        if i != j - 1 and i != j + 1:  # Avoid counting sequential neighbors
                            if self.sequence[j] == 'H':  # Check if the neighbor is also H
                                energy += ENERGY_HH
        return energy // 2  # Each H-H contact is double-counted

                    

    def fold(self, max_iter, temperature,cooling_rate,rules='all', output_energies=False):
        """
        Perform simulated annealing to fold the protein.

        Parameters:
            max_iter (int): Maximum number of iterations.
            temperature (float): Initial temperature for annealing.
            cooling_rate (float): Cooling rate to decrease temperature.
            rule (str): The rule for generating moves ('small', 'rotation', or 'both').

        Returns:
            list of np.ndarray: A list of protein structures representing the folding process.
        """
        if output_energies:
            energies=[]
        energy=self.get_energy()
        states=[self.structure.copy()]
        for i in range(max_iter):
            for move in np.random.permutation(self.get_moves(rules=rules, position=np.random.choice(np.arange(self.structure.shape[0])))):
                new_energy=self.get_energy(move)
                if energy >= new_energy or np.random.rand() < np.exp((energy - new_energy) / max(temperature, 1e-8)):
                    self.structure=move.copy()
                    states.append(move.copy())
                    energy=new_energy
                    if output_energies:
                        energies.append(energy)
                    break
            temperature*=cooling_rate
        if output_energies:
            return states, energies
        return states
    
                    
                    

# Visualisations
def visualize_lattice(positions,sequence,ax=None):
    """
    Visualize the lattice and the protein fold.

    Parameters:
        positions (list of tuple): List of positions of the protein on the lattice.
        sequence (str): The protein sequence.
        ax (matplotlib.axes.Axes, optional): Axes to plot on. If None, a new figure and axes are created.

    Returns:
        tuple: A tuple containing the figure and axes objects.
    """
    if ax is None:
        fig,ax=plt.subplots(1,1)
    x_coords = [pos[0] for pos in positions]
    y_coords = [pos[1] for pos in positions]
    #plt.figure(figsize=(8, 8))
    ax.plot(x_coords, y_coords, marker='o', linestyle='-', color='blue')
    for i, pos in enumerate(positions):
        ax.text(pos[0], pos[1], sequence[i], color="red", fontsize=12, ha='center', va='center')
    ax.set_title("Protein Folding on a 2D Lattice")
    ax.grid(True)
    #ax.set_xlim(-15, 15)  # Fix x-axis range
    #ax.set_ylim(-15, 15)
    return fig, ax


def folding_gif(folding_steps, sequence, output_file="protein_folding.gif", grid_size=None, cell_colors=None, interval=200):
    """
    Save the protein folding process as a GIF.

    Parameters:
        folding_steps (list of np.ndarray): List of intermediate folding states.
        sequence (list of str): The protein sequence.
        output_file (str): Filepath to save the GIF.
        grid_size (tuple, optional): Size of the grid. Defaults to None.
        cell_colors (list, optional): Colors for the cells. Defaults to None.
        interval (int): Time between frames in milliseconds.

    Returns:
        None
    """
    # Initialize the figure
    fig, ax = plt.subplots()
    ax.set_title("Protein Folding on a 2D Lattice")
    ax.grid(True)

    # Create animation frames
    frames = []
    for step in folding_steps:
        elements = []

        # Plot the connections and points
        line, = ax.plot(step[:, 0], step[:, 1], marker='o', linestyle='-', color='blue', animated=True)
        elements.append(line)

        # Annotate each point with its corresponding label
        texts = []
        for i, pos in enumerate(step):
            text = ax.text(pos[0], pos[1], sequence[i], color="red", fontsize=12, ha='center', va='center', animated=True)
            texts.append(text)
        elements.extend(texts)

        # Append elements to the frame
        frames.append(elements)

    # Set axis limits
    all_x = np.concatenate([step[:, 0] for step in folding_steps])
    all_y = np.concatenate([step[:, 1] for step in folding_steps])
    ax.set_xlim(all_x.min() - 1, all_x.max() + 1)
    ax.set_ylim(all_y.min() - 1, all_y.max() + 1)

    # Create the animation
    anim = ArtistAnimation(fig, frames, interval=interval, blit=True)

    # Save the animation as a GIF
    anim.save(output_file, writer="pillow")




import re

def expand_notation(sequence: str) -> str:
    """This is only used to process a format of sequence input from a benchmark dataset

    Args:
        sequence (str): for example: (HP)2PH2PHP2HPH2P2HPH

    Returns:
        str: sequence containing only H and Ps
    """
    while '(' in sequence:
        sequence = re.sub(r'\(([^()]+)\)(\d*)', lambda m: m.group(1) * int(m.group(2) or 1), sequence)
    sequence = re.sub(r'H(\d+)', lambda m: 'H' * int(m.group(1)), sequence)
    sequence = re.sub(r'P(\d+)', lambda m: 'P' * int(m.group(1)), sequence)
    return sequence



def generate_connected_points(n, start=(0, 0)):
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
    
    return sequence
