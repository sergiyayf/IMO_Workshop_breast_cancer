from load.pyMCDS import pyMCDS
from load.pyMCDS_timeseries import pyMCDS_timeseries
import numpy as np
import matplotlib as mpl
#mpl.use('TkAgg')
import matplotlib.pyplot as plt
 
## load data


## Set our z plane and get our substrate values along it
def plot_one_cell_type(cell_df, ax=None, c_type=0, color='grey', edgecolor='dimgrey'):
    x_pos = cell_df['position_x']
    y_pos = cell_df['position_y']
    cell_type = cell_df['cell_type']
    rad = (cell_df['total_volume'] * 3 / 4 / np.pi) ** (1 / 3)
    circles = [plt.Circle((xi, yi), radius=ri, linewidth=0) for xi, yi, ri in
               zip(x_pos[cell_type == c_type], y_pos[cell_type == c_type], rad[cell_type == c_type])]
    c = mpl.collections.PatchCollection(circles, color=color, alpha=1, edgecolor=edgecolor, linewidth=0.01)
    ax.scatter(x_pos[cell_type == c_type], y_pos[cell_type == c_type], color=color, s=1.5, alpha=1, linewidth=0.00,
               edgecolor=edgecolor)

    ax.add_collection(c)
    return
def plot_all_cell_types(ax,mcds):
    cell_df = mcds.get_cell_df()
    cell_types = np.array(cell_df['cell_type'])

    for idx,c_type in enumerate(np.unique(cell_types)):
        plot_one_cell_type(cell_df,ax,c_type,color=plt.get_cmap("tab10").colors[idx])

    xmin = mcds.get_2D_mesh()[0][0][0] - mcds.get_mesh_spacing() / 2
    xmax = mcds.get_2D_mesh()[0][0][-1] + mcds.get_mesh_spacing() / 2
    ymin = mcds.get_2D_mesh()[1][0][0] - mcds.get_mesh_spacing() / 2
    ymax = mcds.get_2D_mesh()[1][-1][-1] + mcds.get_mesh_spacing() / 2

    ax.set_xlim([xmin, xmax])
    ax.set_ylim([ymin, ymax])
    ax.set_aspect('equal')

    return

def plot_microenv(fig,ax=None,cmap=plt.get_cmap('viridis'),vmin=0,vmax=38):

    z_val = 0.00
    chemical = 'drug'; 
    plane_drug = mcds.get_concentrations('drug', z_slice=z_val)

    ## Get the 2D mesh for contour plotting
    xx, yy = mcds.get_2D_mesh()
    # We want to be able to control the number of contour levels so we
    # need to do a little set up
    num_levels = 100

    
    # set up the figure area and add data layers
    #fig, ax = plt.subplots()
    if chemical == 'drug': 
        #min_conc = plane_drug.min()
        #print(min_conc)
        #max_conc = plane_drug.max()
        #my_levels = np.linspace(min_conc, max_conc, num_levels)
        cs = ax.contourf(xx, yy, plane_drug, cmap =cmap,vmin=vmin,vmax=vmax)
        #ax.contour(xx, yy, plane_drug, color='black', levels = my_levels,linewidths=0.5)
        ax.set_title('drug (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(mcds.get_time(),mcds.data['metadata']['time_units'],z_val,mcds.data['metadata']['spatial_units']) )
    elif chemical == 'oxygen': 
        min_conc = plane_oxy.min()
        max_conc = plane_oxy.max()
        my_levels = np.linspace(min_conc, max_conc, num_levels)
        cs = ax.contourf(xx, yy, plane_oxy, levels=my_levels)
        ax.contour(xx, yy, plane_oxy, color='black', levels = my_levels,linewidths=0.5)
        ax.set_title('oxygen (mmHg) at t = {:.1f} {:s}, z = {:.2f} {:s}'.format(mcds.get_time(),mcds.data['metadata']['time_units'],z_val,mcds.data['metadata']['spatial_units']) )
        
    # Now we need to add our color bar
    cbar1 = fig.colorbar(cs, shrink=0.75)
    cbar1.set_label('mmHg')
    
    # Let's put the time in to make these look nice
    ax.set_aspect('equal')
    ax.set_xlabel('x (micron)')
    ax.set_ylabel('y (micron)')


#mcds = pyMCDS('output00000011.xml', '..\\results\\output')
#plot_microenv();
#print(mcds.get_concentrations('drug', z_slice=0.00))

for i in range(0,10):
    mcds = pyMCDS('output0000000'+str(i)+'.xml', r'..\\..\\temp\\PhysiCell_V_1.10.4_0\\output')

    fig, ax = plt.subplots()
    plot_microenv(fig,ax=ax,cmap=plt.get_cmap('inferno'),vmin = 0, vmax = 38)
    #conctr = mcds.get_concentrations('drug', z_slice=0.00)
    #ax.imshow(conctr,cmap = plt.get_cmap('inferno'),vmin = 0, vmax = 38)
    plot_all_cell_types(ax,mcds)


plt.show()

