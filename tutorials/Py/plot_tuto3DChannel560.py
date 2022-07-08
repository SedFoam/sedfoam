from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import os


case = '3DChannel560'
basepath = '../LES/'
sol = basepath + case + '/'

# plot creation and parameters
matplotlib.rcParams.update({'font.size': 20})
matplotlib.rc('text', usetex=True)
plt.figure(figsize=[15, 6])
gs = plt.GridSpec(1, 3)
gs.update(left=0.08, right=0.975, top=0.95, bottom=0.15, wspace=0.1,
          hspace=0.1)
font_legend = matplotlib.font_manager.FontProperties(family='Comic Sans MS',
                                                     weight='bold',
                                                     style='normal', size=14)

h = 0.02

# plot 1
ax1 = plt.subplot(gs[0, 0])
ax1.axis([1e-10, 0.3, 0, 1])
ax1.set_xlabel(r'$\langle\bar\phi\rangle$')
ax1.set_ylabel(r'$y/h$')
ax1.grid(linestyle=':')

# plot 2
ax2 = plt.subplot(gs[0, 1])
ax2.axis([0, 0.65, 0, 1])
ax2.set_xlabel(r'$\langle\tilde u^f\rangle_F$, $\langle\tilde u^s\rangle_F$')
ax2.set_yticklabels([])
ax2.grid(linestyle=':')

# plot 3
ax3 = plt.subplot(gs[0, 2])
ax3.axis([0, 1e-2, 0, 1])
ax3.set_xlabel(r"$\langle\tilde u^{f'} \tilde v^{f'}\rangle_F$, "
               r"$\langle\tilde u^{s'} \tilde v^{s'}\rangle_F$")
ax3.set_yticklabels([])
ax3.grid(linestyle=':')

# Numerical data reading
avg_file = Dataset(os.path.join(sol, 'postProcessing', case+'_averaged.nc'))
phi = avg_file.variables['alpha_a'][:]
umean_s = avg_file.variables['uaf_a'][:]
umean_f = avg_file.variables['ubf_a'][:]
uprim_s = avg_file.variables['uaprimf_a'][:]
uprim_f = avg_file.variables['ubprimf_a'][:]
y = avg_file.variables['y'][:]

# Plot data
ax1.semilogx(phi, y/h, '-k')
ax2.plot(umean_s[0], y/h, '-r', label='solid')
ax2.plot(umean_f[0], y/h, '-k', label='fluid')
ax3.plot(uprim_s[0], y/h, '-r', label='solid')
ax3.plot(uprim_f[0], y/h, '-k', label='fluid')

ax2.legend()
ax3.legend()

plt.savefig('Figures/tuto3DChannel560.png', facecolor='w', edgecolor='w',
            format='png')
plt.show()
