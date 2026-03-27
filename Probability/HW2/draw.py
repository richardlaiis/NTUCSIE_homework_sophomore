import matplotlib.pyplot as plt
import numpy as np

# 1. PMF Data for X = min(d1, d2)
# f(x) = [2(7-x)-1] / 36
x_values = np.arange(1, 7)
x_probs = [(2 * (7 - x) - 1) / 36 for x in x_values]

# 2. PMF Data for Y = |d1 - d2|
# g(y) is 6/36 for y=0, and 2(6-y)/36 for y=1,2,3,4,5
y_values = np.arange(0, 6)
y_probs = [6/36 if y == 0 else (2 * (6 - y)) / 36 for y in y_values]

def save_pmf_histogram(values, probs, title, filename, color):
    plt.figure(figsize=(8, 6))
    
    # width=1.0 makes the bars touch, creating a proper histogram look
    plt.bar(values, probs, color=color, edgecolor='black', width=1.0, alpha=0.8)
    
    plt.title(title, fontsize=14)
    plt.xlabel('Value', fontsize=12)
    plt.ylabel('Probability P(Value)', fontsize=12)
    plt.xticks(values)
    
    # Adding text labels on top of bars for clarity
    for i, p in enumerate(probs):
        plt.text(values[i], p + 0.005, f'{p:.3f}', ha='center', va='bottom')

    plt.grid(axis='y', linestyle='--', alpha=0.3)
    
    # Save the file
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close() # Close the figure to free up memory
    print(f"Successfully saved: {filename}")

# Execute saves
save_pmf_histogram(x_values, x_probs, 'PMF of X (Minimum of Two Dice)', 'hist_X_min.png', 'skyblue')
save_pmf_histogram(y_values, y_probs, 'PMF of Y (Range of Two Dice)', 'hist_Y_range.png', 'salmon')