from matplotlib.patches import Circle
from casa_cube import casa_cube as casa
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pymcfost as mcfost
from my_utils import colorbar2, shift_axes, get_channel_maps, colorbar

fontsize = 12

plt.rcParams['axes.labelsize']  = fontsize
plt.rcParams['legend.fontsize'] = fontsize
plt.rcParams["figure.figsize"]  = [10.,10.]
plt.rcParams['xtick.labelsize'] = fontsize
plt.rcParams['ytick.labelsize'] = fontsize
plt.rcParams['axes.titlesize']  = fontsize
plt.rcParams['lines.dash_capstyle'] = 'round'
plt.rcParams['lines.linewidth'] = 3

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


#------------------------------#
# Directory Information
#------------------------------#
# Simulation name
name = "kinks2_fig1_taper_outerrad"
run = "HD_163296_kinks2"

# Path to the observations data
dir = f"/home/tom/Documents/Protoplanetary_Discs_Project/wakeflow/{run}/{name}/"

# Model directories
mod_dir_basename = dir

# Name of the output PDF file
filename = "{}.pdf".format(name)
#------------------------------#




#------------------------------#
# Observational Variables
#------------------------------#

# Include Observation Data
include_observation = True

# Directory to Observational Data
# Velocity CO channel must be called lines.fits.gz or lines.fits
# Continuum data must be called RT.fits.gz or RT.fits
#dir_observation = dir + "Observation/"
dir_observation = "/home/tom/Documents/Protoplanetary_Discs_Project/spiralmoment/line/maps/lines/"

# The name of the observation
obs_name = "HD 163296"

# Match Observation Scales
# This will match the bmin, bmax and bpa from the observation
match_observation = True
#------------------------------#



#------------------------------#
# System Variables
#------------------------------#

# Star velocity factor moving away
v_system = 5.855

# Location of planet
# For simulations without a planet, use 0, 0
p_loc = [0, 0]

# Location of star shift
s_shift = [0, 0]

# System distance
distance = 101.5  #  ---  only used for 2D analytics plots
#------------------------------#



#------------------------------#
# Model Variables
#------------------------------#

# The model names

models = ["Simulation", "5.0Mj", "3.0Mj", "1.0Mj"]
model_labels = [
    r"2 $\mathrm{M_J}$ Sim", r"5 $\mathrm{M_J}$", r"3 $\mathrm{M_J}$", r"1 $\mathrm{M_J}$"
]

models = ["Simulation", "5.0Mj"]
model_labels = [
    r"Simulation", r"Analytics"
]

"""
models = ["5.0Mj"]
model_labels = [
    r"A: 5 $\mathrm{M_J}$"
]
"""

"""
models = ["3Mj"]
model_labels = [
    r"A: 3 $\mathrm{M_J}$"
]
"""

# The mass of the planets (in Jupiter masses)
# For models without planets, write 0
# MAKE SURE THE SIZE OF THIS ARRAY IS THE SAME AS MODELS
p_masses = [2,5,3,1]
p_masses = [2,5]
#p_masses = [5]

# Velocity channels
#v_channels = [0.84, 0.96, 1.08, -0.96]
v_channels = [-0.4, -0.2, 0.2, 0.4]
#v_channels = [0.4, 0.7, 1, 1.3, 1.6]

# Include 2D analytics?
include_2D_analytics = False

# 2D analytics name
name_2D = "2D"

# 2D analytics label
label_2D = "2D A"

# 2D analytics channel resolution [km/s]
ch_res_2d = 0.05
#------------------------------#



#------------------------------#
# Plotting Variables
#------------------------------#

# label coords
label_x_coord = 0.1
label_y_coord = 0.82

# Whether to plot continuum images or not
include_continuum = False

# The colour mapping for the observational plot
cmap_cont = "gist_earth"

# Whether to use Flux Temperature or not
c_plot_temp = True

# Whether to plot the location of sinks on the data or not
plot_sinks = False

# The minimum and maximum flux value for the pixels for velocity images
# These will scale the velocity channels
v_f_min = 9
v_f_max = 80

# The f_min value for the pixels in the continuum image
# The minimum and maximum flux value for the pixels on the continuum images
# These will scale the continuum
c_f_min = 2
c_f_max = 50

# The continuum pixel addition
# This makes the continuum darker the larger the value is
# Keep this as 0 if you do not want to change the scaling
c_mod_pix_add = 0

# The continuum colour scale
# Options are 'log' or 'lin'
c_color_scale = 'lin'

# The limits of the graph
# (max_x, min_x, min_y, max_y)
lim_num = 5
limits = [lim_num, -lim_num, -lim_num, lim_num]
#limits = [4.5,-3,-3, 4.5]
#limits = [500,-500,-500,500]

x_ticks = [-2.5, 0, 2.5]
y_ticks = [-2.5, 0, 2.5]

# Python Figure Size
f_size = 2.0

# Python spacing between figures
f_spacing = 0.0
#------------------------------#




#------------------------------#
# Applicaion Variables
# DO NOT CHANGE
#------------------------------#

# The number of models being used
n_models = len(models)

# The number of channels being used
n_channels = len(v_channels)
#------------------------------#




#------------------------------#
# Set up the Plots
#------------------------------#

# If including the observational data
if include_observation:
    # Create continuum and CO data
    if include_continuum:
        cont =  casa.Cube(dir_observation + "RT.fits")

    CO =        casa.Cube(dir_observation + "lines.fits")

# Create the subplots
fig, axes = plt.subplots(
    nrows = n_models + int(include_observation) + int(include_2D_analytics),
    ncols = n_channels + int(include_continuum),
    figsize = (f_size * (n_channels + int(include_continuum)), f_size * (n_models + int(include_observation) + int(include_2D_analytics))),
    sharex='all',
    sharey='all'
)

for i in range(int(include_observation)+len(models)):
    for j in range(n_channels):
        axes[i,j].set_xticks(x_ticks)
        axes[i,j].set_yticks(y_ticks)
        axes[i,j].xaxis.set_ticks_position('both')
        axes[i,j].yaxis.set_ticks_position('both')
        axes[i,j].xaxis.set_tick_params(direction='in', color='white', labelcolor='black')
        axes[i,j].yaxis.set_tick_params(direction='in', color='white', labelcolor='black')
        axes[i,j].spines['bottom'].set_color('white')
        axes[i,j].spines['top'].set_color('white') 
        axes[i,j].spines['right'].set_color('white')
        axes[i,j].spines['left'].set_color('white')

# Add some whitespace between them
plt.subplots_adjust(wspace = f_spacing, hspace = f_spacing)

# Adjust the axes for each row
if include_continuum:
    for i in range(n_models + int(include_observation)):
        shift_axes(axes[i,0],-0.03,0)
        #shift_axes(axes[i,4:],0.01,0)

#------------------------------#




#------------------------------#
# Creates a circle at the planet position on the graph
def CreateCircle ():
    return Circle(
        (p_loc[0], p_loc[1]),
        0.15,
        clip_on = False,
        zorder = 10,
        linewidth = 2,
        edgecolor = 'white',
        linestyle = ":",
        facecolor = (0, 0, 0, .0125),
        alpha = 0.3,
        fill = False
    )
#------------------------------#




#------------------------------#
# Broadcast any Errors
#------------------------------#

# Check for mismatched arrays
if len(p_masses) != len(models):
    raise Exception("Incorect Planet Masses Array Size. Make sure the array length is identical to the model array length.")




#------------------------------#
# Plot the Observational data
#------------------------------#

# If including the observational data
if include_observation:

    # If plotting continuum graphs
    if include_continuum:

        print("Plotting Observation Continuum")

        # We plot the observations on the first row
        image = cont.plot(
            colorbar = False,
            cmap = cmap_cont,
            color_scale = c_color_scale,
            ax = axes[0,0],
            no_xlabel = True,
            no_ylabel = False,
            limits = limits,
            shift_dx = s_shift[0],
            shift_dy = s_shift[1],
            Tb = c_plot_temp,
            fmin = c_f_min,
            fmax = c_f_max
        )

        # Show the colour bar in the plot
        #colorbar2(image)

        # Add the planet and star to the plot
        axes[0,0].plot(s_shift[0],  s_shift[1],     "*", color="white", ms=4)
        axes[0,0].plot(p_loc[0],    p_loc[1],       "o", color="cyan",  ms=2)

        # Label the planet name on the continuum plot
        axes[0,0].text(
            label_x_coord,
            label_y_coord,
            obs_name,
            horizontalalignment = 'left',
            color = "white",
            transform = axes[0,0].transAxes,
            fontsize = 10
        )


    print("Plotting Observation Velocity Channels")

    # Loop though all the channels
    for i in range(n_channels):

        # Determine the current velocity for the observational data
        iv = np.abs(CO.velocity - (v_system + v_channels[i])).argmin()

        # Only show color bar in the last channel
        show_colorbar = i == n_channels - 1

        # Only show y-label for first plot
        if i == 0:
            dont_show_y_label = False
        else:
            dont_show_y_label = True

        # Plot the velocity channel
        vel_im = CO.plot(
            iv = iv,  
            v0 = v_system,
            colorbar = False,
            ax = axes[0, int(include_continuum) + i],
            no_xlabel = True,
            no_ylabel = dont_show_y_label,
            limits = limits,
            shift_dx = s_shift[0],
            shift_dy = s_shift[1],
            Tb = c_plot_temp,
            fmax = v_f_max,
            fmin = v_f_min
        )

        if show_colorbar:
            #colorbar2(vel_im,  label=r"$T_b$ [K]")
            pass
        
        axes[0, 0].text(
            label_x_coord,
            label_y_coord,
            "MAPS",
            horizontalalignment = 'left',
            color = "white",
            transform = axes[0,0].transAxes,
            fontsize = 10
        )

        # Add a circle where the planet is expected to be
        circle = CreateCircle()
        axes[0, int(include_continuum)  + i].add_artist(circle)
#------------------------------#




#------------------------------#
# Plot the 2D Analytic Model
#------------------------------#

if include_2D_analytics:

    # get grid for analytics in au
    r = np.linspace(1,700,1000)
    phi = np.linspace(0,2*np.pi,1000)
    R,PHI = np.meshgrid(r,phi)
    X = R*np.cos(PHI)
    Y = R*np.sin(PHI)
    distance *= 206265     # convert pc to au

    # convert to angular coordinates using distance
    X_ang = 206265 * np.arctan(X / distance)
    Y_ang = 206265 * np.arctan(Y / distance)

    # get vz_field from file
    full_vz_field = np.load(f"{dir}{name_2D}/vz_field.npy")
    vz_field = full_vz_field[:,:,2]

    # plot each channel
    for i in range(n_channels):

        channel = [v_channels[i]]

        # get channel maps
        channel_maps = get_channel_maps(vz_field, channel, ch_res_2d)

        # set channel maps levels
        levelsch = [channel[0] - ch_res_2d]
        Nvch = len(channel)
        levelsch += [0.5*(channel[j]+channel[j+1]) for j in range(0,Nvch-1)]
        levelsch += [channel[-1] + ch_res_2d]
        levelsch = np.array(levelsch)

        # plot channel maps
        chmp = axes[0 + int(include_observation), i + int(include_continuum)].contourf(X, Y, channel_maps, levels = levelsch, colors = 'orange', zorder = 2)
        cm = plt.contour(X, Y, channel_maps, levels = levelsch, colors='black', linewidths = 1 ,zorder = 3)

        # legend
        proxy = [plt.Rectangle((0,0),1,1,fc=pc.get_facecolor()[0]) for pc in chmp.collections]
        leg2 = axes[0 + int(include_observation), i + int(include_continuum)].legend(proxy, [str(v_channels[j]) for j in range(len(v_channels))],loc="upper left",title="$v_{ch}$ [km/s]:")

        # axis label
        if i == 0:
            axes[0 + int(include_observation), i + int(include_continuum)].set_ylabel(r'$\Delta$ Dec ["]')

    axes[0 + int(include_observation), 0].text(
        label_x_coord,
        label_y_coord,
        label_2D,
        horizontalalignment = 'left',
        color = "black",
        fontsize = 10
    )
#------------------------------#




#------------------------------#
# Plot the Simulation Models
#------------------------------#

# Loop through each of the simulaions
for k, mod in enumerate(models):
    # Get the model directory
    mod_dir = mod_dir_basename + str(mod) + "/data_CO"
    
    print(mod_dir)

    # Print message
    print("Analysing output from {} located at:\n\t{}".format(mod, mod_dir))

    # Get the continuum and CO images
    #mod_cont = mcfost.Image(mod_dir)
    mod_CO = mcfost.Line(mod_dir)

    # Determine whether to display the xlabel or not
    no_xlabel = k < n_models - 1

    # Add pixel values to the continuum image
    #mod_cont.image += c_mod_pix_add * mod_cont.image[4,0,0,:,:]

    # If plotting the continuum
    if include_continuum:

        # Plot the continuum
        if include_observation and match_observation:
            image = mod_cont.plot(
                ax = axes[k + 1, 0],
                colorbar = False,
                bmaj = cont.bmaj,
                bmin = cont.bmin,
                bpa = cont.bpa,
                no_xlabel = no_xlabel,
                limits = limits,
                cmap = cmap_cont,
                scale = c_color_scale,
                Tb = c_plot_temp,
                vmin = c_f_min,
                vmax = c_f_max,
                plot_stars = plot_sinks
            )
        
        # If no observational data to base scales on
        else:
            image = mod_cont.plot(
                ax = axes[k + int(include_observation), 0],
                colorbar = False,
                no_xlabel = no_xlabel,
                limits = limits,
                cmap = cmap_cont,
                scale = c_color_scale,
                Tb = c_plot_temp,
                vmin = c_f_min,
                vmax = c_f_max,
                plot_stars = plot_sinks
            )

        # Label the model on the continuum plot
        if p_masses[k] == 0:
            mod_label = models[k] + " Model"
        else:
            mod_label = str(p_masses[k]) + "M$_\mathrm{jup}$ Model"
        
        axes[k + int(include_observation), 0].text(
            label_x_coord,
            label_y_coord,
            mod_label,
            horizontalalignment = 'left',
            color = "white",
            transform = axes[k+int(include_observation)+int(include_2D_analytics),0].transAxes,
            fontsize = 10
        )

        # Add the planet and star to the plot
        #axes[k+int(include_observation),0].plot(s_shift[0],  s_shift[1], "*", color="white", ms=4)

        #if mplanet != 0:  
        #    axes[k+int(include_observation),0].plot(p_loc[0],    p_loc[1],   "o", color="cyan",  ms=2)

        # Show the colour bar in the plot
        #colorbar2(image)

    # Set the change in velocity
    delta_v = None

    # Loop through each of the channels
    for i in range(n_channels):

        # Only show color bar in the last channel
        show_colorbar = i == n_channels - 1

        # Add y-axis label for first plot
        if i==0:
            plot_y_label = False
        else:
            plot_y_label = True

        # Plot the CO velocity channel
        if include_observation and match_observation:
            vel_im = mod_CO.plot_map(
                v = v_channels[i], 
                ax = axes[k + 1 + int(include_2D_analytics), int(include_continuum)  + i],
                #ax = axes[int(include_continuum)  + i],
                colorbar = False,
                bmaj = CO.bmaj,
                bmin = CO.bmin,
                bpa = CO.bpa,
                no_xlabel = no_xlabel,
                no_ylabel = plot_y_label,
                limits = limits,
                Tb = c_plot_temp,
                Delta_v = delta_v,
                fmax = v_f_max,
                fmin = v_f_min,
                plot_stars = plot_sinks
            )

        # If no observational data to base scales on
        else:
            vel_im = mod_CO.plot_map(
                v = v_channels[i], 
                ax = axes[k + int(include_observation) + int(include_2D_analytics), int(include_continuum) + i],
                #ax = axes[int(include_continuum) + i],
                colorbar = False,
                no_xlabel = no_xlabel,
                no_ylabel = plot_y_label,
                limits = limits,
                Tb = c_plot_temp,
                Delta_v = delta_v,
                fmax = v_f_max,
                fmin = v_f_min,
                plot_stars = plot_sinks
            )
        
        mod_label = model_labels[k]

        # Plot the circle where the planet is expected to be
        if p_masses[k] != 0:
            circle = CreateCircle()
            axes[k + int(include_observation) + int(include_2D_analytics), int(include_continuum) + i].add_artist(circle)
            #axes[int(include_continuum) + i].add_artist(circle)

        axes[k + int(include_observation) + int(include_2D_analytics), 0].text(
        #axes[0].text(
            label_x_coord,
            label_y_coord,
            mod_label,
            horizontalalignment = 'left',
            color = "white",
            transform = axes[k+int(include_observation)+int(include_2D_analytics),0].transAxes,
            #transform = axes[0].transAxes,
            fontsize = 10
        )

        # Add in a colorbar
        if show_colorbar:
            #colorbar2(vel_im, label=r"$T_b$ [K]")
            pass

fig.subplots_adjust(right=0.9)
cbar_ax = fig.add_axes([0.902, 0.111, 0.02, 0.768])
cb = fig.colorbar(vel_im, cax=cbar_ax, label=r"$T_b$ [K]")

cbar_ax.yaxis.set_ticks_position('none')
#cbar_ax.yaxis.set_tick_params(direction='in', color='black', labelcolor='black')

#cb.outline.set_edgecolor('black')

# plot dotted circles
for i in range(3):
    
    size = 220
    alpha = 0.3
    
    # first column
    axes[i,0].scatter(2.35, 1.9, s=size+40,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    #axes[i,0].scatter(0, -1.25, s=size,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    
    # second
    axes[i,1].scatter(1, 1.5, s=size+20,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    #axes[i,1].scatter(0, -1.25, s=size,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    
    # third
    axes[i,2].scatter(0.2, 2.05, s=size,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    axes[i,2].scatter(-2.6, -2.5, s=size,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    
    # fourth
    axes[i,3].scatter(0.2, 2.05, s=size,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)
    axes[i,3].scatter(-2.75, -2.2, s=size+20,  facecolors='none', edgecolors='white', linestyle='--', alpha=alpha)


#------------------------------#

for i in range(int(include_observation)+len(models)):
    for j in range(n_channels):
        axes[i,j].set_xticks(x_ticks)
        axes[i,j].set_yticks(y_ticks)

# Save the figure
plt.savefig("secondary_" + filename, bbox_inches='tight')

# Show the graph in an xw display window
plt.show()
