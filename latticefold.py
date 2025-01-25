import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation

# Define constants
ENERGY_HH = -1  # Energy for H-H contact
GRID_SIZE = 20  # Size of the lattice grid


# Utils
rotation_matrix = lambda theta: (np.array([[np.cos(theta), -np.sin(theta)],[np.sin(theta), np.cos(theta)]]).astype(int))
def rotation(shape,center=np.array([0,0]).reshape((2,1)),theta=0):
    return np.dot(rotation_matrix(theta),shape-center.reshape((2,1)))+center.reshape((2,1))

def are_neighbors(p1, p2):
    """
    Check if two points are neighbors on a 2D lattice.
    """
    x1, y1 = p1
    x2, y2 = p2
    return abs(x1 - x2) + abs(y1 - y2) == 1

def get_neighbors(position):
    """Get the coordinates of all neighbors of a given position."""
    x, y = position
    return [(x + 1, y), (x - 1, y), (x, y + 1), (x, y - 1)]

def sterical_clash(pos):
    """Checks for sterical clash in a position"""
    return not len(np.unique(pos, axis=0))==len(pos)

def chain_connected(strucutre):
    """
    Checks if a chain of points in a 2D lattice is connected.

    Args:
        points (list of tuple): List of points representing the chain, where each point is a tuple (x, y).

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
            self.structure=np.c_[np.ones(len(sequence),dtype=int),np.arange(len(sequence))]
        else: self.structure=initial_structure


    def get_moves(self, check=True, rule='small'):
        moves=[]
        if rule == 'rotation' or rule == 'both':
            for i,pos in enumerate(self.structure):
                if i>0:
                    for theta in np.arange(1,4)/2*np.pi:
                        move=self.structure.T.copy()
                        move[:,i:]=rotation(move[:,i:], pos.reshape(2,1),theta)
                        if not check:
                            moves.append(move.T)
                        elif chain_connected(move.T) and not sterical_clash(move.T):
                            moves.append(move.T)
                    
        if rule == 'small' or rule == 'both':
            for nb in get_neighbors(self.structure[1]):
                if not (nb==self.structure[1]).all():
                    move=self.structure.copy()
                    move[0]=nb
                    if not check:
                        moves.append(move)
                    elif chain_connected(move) and not sterical_clash(move):
                        moves.append(move.copy())

            for nb in get_neighbors(self.structure[-2]):
                if not (nb==self.structure[-2]).all():
                    move=self.structure.copy()
                    move[-1]=nb
                    if not check:
                        moves.append(move)
                    elif chain_connected(move) and not sterical_clash(move):
                        moves.append(move.copy())

            for i in range(1,len(self.structure)-1):
                move=self.structure.copy()
                if np.abs(move[i-1]-move[i+1])[0]==1:
                    move[i]=move[i-1]-move[i]+move[i+1]
                    if not check:
                        moves.append(move)
                    elif chain_connected(move) and not sterical_clash(move):
                        moves.append(move)
                    
        return moves

    def get_energy(self, structure=None):
        """Calculate the total energy based on H-H contacts."""
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

                    

    def fold(self, max_iter, temperature,cooling_rate,rule):
        energy=self.get_energy()
        states=[self.structure.copy()]
        for i in range(max_iter):
            for move in np.random.permutation(self.get_moves(rule=rule)):
                new_energy=self.get_energy(move)
                if energy >= new_energy or np.random.rand() > np.exp(energy - new_energy) / max(temperature, 1e-8):
                    self.structure=move.copy()
                    states.append(move.copy())   # this may need a .copy()
                    energy=new_energy
                    break
            temperature*=cooling_rate
        return states
    
                    
                    

# Visualisations
def visualize_lattice(positions,sequence,ax=None):
    """Visualize the lattice and the protein fold."""
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
    Saves the lattice protein folding steps as a GIF using matplotlib's ArtistAnimation.

    Parameters:
        folding_steps (list of np.ndarray): List of 2D numpy arrays representing intermediate folding states.
        sequence (list of str): Sequence of labels corresponding to each point in the folding steps.
        output_file (str): Filepath to save the GIF.
        interval (int): Time between frames in milliseconds.
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
