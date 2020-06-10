# linear regression and correlations
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt


x = np.arange(10, 20)
y = np.array([2, 1, 4, 5, 8, 12, 18, 25, 96, 48])

# Pearson's correlation coefficient (r)
r, ppval = scipy.stats.pearsonr(x, y)
print(r, ppval, sep='\t')

# Spearman's correlation coefficient (ro)
ro, spval = scipy.stats.spearmanr(x, y)
print(ro, spval, sep='\t')

# Kendall's correlation coefficient (tau)
tau, kpval = scipy.stats.kendalltau(x, y)
print(tau, kpval, sep='\t')

# Linear regression and Pearson's r
slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)

# plot
plt.style.use('ggplot')
line = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r={r:.2f}'
fig, ax = plt.subplots()
ax.plot(x, y, linewidth=0, marker='o', label='Data points')
ax.plot(x, intercept + slope * x, label=line)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.legend(facecolor='white')
plt.show()
