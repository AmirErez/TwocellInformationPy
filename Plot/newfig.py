import matplotlib.pyplot as plt
from matplotlib import rc



def newfig(width_cm=4.3, height_cm=2.65, font_size=10, nrows=1, ncols=1):
    """
    Create a publication-quality figure using matplotlib.

    Parameters:
    - width_cm: float, width of the figure in centimeters
    - height_cm: float, height of the figure in centimeters (optional)
    - font_size: int, font size for text elements

    Returns:
    - fig: matplotlib.figure.Figure, figure handle
    - axs: array of matplotlib.axes._subplots.AxesSubplot, axis handle
    """

    # Enable LaTeX interpretation and set font family
    rc('text', usetex=True)
    plt.rc('font', family='serif')

    # Convert width and height from cm to inches
    width_in = width_cm / 2.54
    height_in = height_cm / 2.54

    # Create figure and axis handles with high-resolution output
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(width_in, height_in), dpi=300)
    # Set font size and line width for better readability
    plt.rc('font', size=font_size)
    plt.rc('lines', linewidth=1.5)
    return fig, axs