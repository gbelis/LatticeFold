import matplotlib.pyplot as plt
from matplotlib.animation import ArtistAnimation
import numpy as np

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