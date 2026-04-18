import numpy as np
import matplotlib.pyplot as plt

# ── Helper: draw a clean panel with vertical/horizontal guides ──────────────
def style_ax(ax, title):
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.axhline(0, color='k', linewidth=0.7)   # horizontal zero line
    ax.axvline(0, color='k', linewidth=0.7)   # vertical zero line
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=9)

fig, axes = plt.subplots(3, 2, figsize=(10, 12))
fig.suptitle("Problem 3.1-8: PDF and CDF Sketches", fontsize=14, fontweight='bold')

# ══════════════════════════════════════════════════════
# PART (a): f(x) = x³/4,  F(x) = x⁴/16,  0 < x < 2
# ══════════════════════════════════════════════════════
x_a = np.linspace(0, 2, 300)      # 300 evenly-spaced points in [0, 2]
f_a = x_a**3 / 4                  # pdf formula
F_a = x_a**4 / 16                 # cdf formula

ax = axes[0, 0]
ax.plot(x_a, f_a, color='steelblue', label='$f(x)=x^3/4$')
ax.fill_between(x_a, f_a, alpha=0.15, color='steelblue')   # shaded area
# mark mu
mu_a = 8/5
ax.axvline(mu_a, color='red', linestyle='--', label=f'$\\mu = 8/5 = {mu_a:.2f}$')
style_ax(ax, '(a) PDF: $f(x) = x^3/4$')
ax.set_xlabel('x'); ax.set_ylabel('f(x)')

ax = axes[0, 1]
ax.plot(x_a, F_a, color='darkorange', label='$F(x)=x^4/16$')
ax.axhline(1, color='gray', linestyle=':')      # asymptote at F=1
style_ax(ax, '(a) CDF: $F(x) = x^4/16$')
ax.set_xlabel('x'); ax.set_ylabel('F(x)')

# ══════════════════════════════════════════════════════
# PART (b): f(x) = (3/16)x²,  F(x) = (x³+8)/16,  -2 < x < 2
# ══════════════════════════════════════════════════════
x_b = np.linspace(-2, 2, 300)
f_b = (3/16) * x_b**2
F_b = (x_b**3 + 8) / 16

ax = axes[1, 0]
ax.plot(x_b, f_b, color='steelblue', label='$f(x)=\\frac{3}{16}x^2$')
ax.fill_between(x_b, f_b, alpha=0.15, color='steelblue')
ax.axvline(0, color='red', linestyle='--', label='$\\mu = 0$ (symmetry)')
style_ax(ax, '(b) PDF: $f(x) = (3/16)x^2$')
ax.set_xlabel('x'); ax.set_ylabel('f(x)')

ax = axes[1, 1]
ax.plot(x_b, F_b, color='darkorange', label='$F(x)=\\frac{x^3+8}{16}$')
ax.axhline(1, color='gray', linestyle=':')
style_ax(ax, '(b) CDF: $F(x) = (x^3+8)/16$')
ax.set_xlabel('x'); ax.set_ylabel('F(x)')

# ══════════════════════════════════════════════════════
# PART (c): f(x) = 1/(2√x),  F(x) = √x,  0 < x < 1
# Note: start from small ε>0 to avoid division by zero at x=0
# ══════════════════════════════════════════════════════
x_c = np.linspace(0.001, 1, 500)   # avoid x=0 (singularity)
f_c = 1 / (2 * np.sqrt(x_c))
F_c = np.sqrt(x_c)

ax = axes[2, 0]
ax.plot(x_c, f_c, color='steelblue', label='$f(x)=1/(2\\sqrt{x})$')
ax.fill_between(x_c, f_c, alpha=0.15, color='steelblue')
ax.axvline(1/3, color='red', linestyle='--', label='$\\mu = 1/3$')
ax.set_ylim(0, 8)                  # clip the spike near x=0 for readability
# annotate the unbounded behavior
ax.annotate('Unbounded\nas $x→0^+$', xy=(0.05, 6),
            xytext=(0.25, 6.5), fontsize=8,
            arrowprops=dict(arrowstyle='->', color='black'))
style_ax(ax, '(c) PDF: $f(x) = 1/(2\\sqrt{x})$')
ax.set_xlabel('x'); ax.set_ylabel('f(x)')

ax = axes[2, 1]
ax.plot(x_c, F_c, color='darkorange', label='$F(x)=\\sqrt{x}$')
ax.axhline(1, color='gray', linestyle=':')
style_ax(ax, '(c) CDF: $F(x) = \\sqrt{x}$')
ax.set_xlabel('x'); ax.set_ylabel('F(x)')

plt.tight_layout()
plt.savefig('hw3_problem4.png', dpi=150, bbox_inches='tight')
plt.show()