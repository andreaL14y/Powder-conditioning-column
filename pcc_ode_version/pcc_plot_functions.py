import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider

def plot_sections_over_time(
        moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, height_of_interest,
        n_space_steps, discrete_time, moisture_gas_initial_bed, moisture_gas_initial_in, moisture_particle_initial,
        moisture_particle_saturated, temp_min, kelvin, hours, max_temp_gas, max_temp_particle):
    epsilon = 0.001
    height_of_interest -= 1
    fig, ax = plt.subplots(n_space_steps, 4, figsize=(20, 13))
    fig.suptitle(f'Moisture & temperature in cylinder sections over time. Total time: {int(hours)} hours.', fontsize=16)
    ax[0, 0].set_title('Moisture Gas')
    ax[0, 1].set_title('Moisture Particle')
    ax[0, 2].set_title('Temp Gas')
    ax[0, 3].set_title('Temp Particle')

    # Set styles
    moisture_color = 'navy'
    temp_color = 'orangered'
    initial_line = 'dashed'
    initial_color = 'gray'
    saturated_color = 'lightcoral'

    for step in range(n_space_steps):
        ax[step, 0].set_ylabel(f'Section {step + 1}                      ', rotation=0, size='large')
        ax[step, 0].plot(discrete_time, moisture_gas_vector[:, height_of_interest, step], c=moisture_color)
        patch = mpatches.Patch(color=moisture_color, label=f'Y {step}')
        ax[step, 0].set_ylim(moisture_gas_initial_bed - epsilon, moisture_gas_initial_in + epsilon)
        ax[step, 0].hlines(moisture_gas_initial_bed, 0, discrete_time[-1], colors=initial_color,
                           linestyles=initial_line)
        ax[step, 0].text(hours * 4 / 5, moisture_gas_initial_bed + epsilon,
                         ('{:.4f}'.format(moisture_gas_initial_bed)), ha='left', va='center')
        ax[step, 0].hlines(moisture_gas_initial_in, 0, discrete_time[-1], colors=saturated_color,
                           linestyles=initial_line)
        ax[step, 0].text(1, moisture_gas_initial_in - epsilon, ('{:.4f}'.format(moisture_gas_initial_in)),
                         ha='left', va='center')


        ax[step, 1].plot(discrete_time, moisture_particle_vector[:, height_of_interest, step],
                         c=moisture_color)
        patch = mpatches.Patch(color=moisture_color, label=f'X {step}')

        # Add lines max and min
        ax[step, 1].hlines(moisture_particle_initial, 0, discrete_time[-1], colors=initial_color,
                           linestyles=initial_line)
        ax[step, 1].text(hours * 4 / 5, moisture_particle_initial + 2 * epsilon,
                         ('{:.4f}'.format(moisture_particle_initial)), ha='left', va='center')

        ax[step, 1].hlines(moisture_particle_saturated, 0, discrete_time[-1], colors=saturated_color,
                           linestyles=initial_line)
        ax[step, 1].text(1, moisture_particle_saturated - 2 * epsilon,
                         ('{:.4f}'.format(moisture_particle_saturated)),
                         ha='left', va='center')


        ax[step, 2].plot(discrete_time, temp_gas_vector[:, height_of_interest, step], c=temp_color)
        patch = mpatches.Patch(color=temp_color, label=f'{step}')
        ax[step, 2].set_ylim(temp_min - 1, max_temp_gas + 1)
        ax[step, 2].hlines(temp_min, 0, discrete_time[-1], colors=initial_color,
                           linestyles=initial_line)
        ax[step, 2].hlines(max_temp_gas, 0, discrete_time[-1], colors=saturated_color,
                           linestyles=initial_line)
        ax[step, 2].text(hours * 4 / 5, max_temp_gas - 0.7,
                         ('{:.2f}'.format(max_temp_gas)), ha='left', va='center')


        ax[step, 3].plot(discrete_time, temp_particle_vector[:, height_of_interest, step],
                         c=temp_color)
        patch = mpatches.Patch(color=temp_color, label=f'{step}')
        ax[step, 3].set_ylim(temp_min - 1, max_temp_particle + 1)
        ax[step, 3].hlines(temp_min, 0, discrete_time[-1], colors=initial_color,
                           linestyles=initial_line)

        ax[step, 3].hlines(max_temp_particle, 0, discrete_time[-1], colors=saturated_color,
                           linestyles=initial_line)
        ax[step, 3].text(hours * 4 / 5, max_temp_particle - 0.7,
                         ('{:.2f}'.format(max_temp_particle)), ha='left', va='center')

        # ax[step, feature].legend(handles=[patch], loc="lower center")
        ax[step, 0].grid()
        ax[step, 1].grid()
        ax[step, 2].grid()
        ax[step, 3].grid()
    plt.savefig('system_over_time.pdf')
    # plt.show()

def plot_heatmap(
        moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, height_of_interest,
        n_space_steps, discrete_time, moisture_gas_initial_bed, moisture_gas_initial_in, moisture_particle_initial,
        moisture_particle_saturated, temp_min, kelvin, hours, max_temp_gas, max_temp_particle):
    time, height, length = np.shape(moisture_gas_vector)
    n_plots = 5
    plot_times = np.linspace(0, time-1, n_plots, dtype=int)
    fig, ax = plt.subplots(n_plots, 4, figsize=(20, 13))

    fig.suptitle(f'Moisture & temperature in cylinder at chosen times. Total time: {int(hours)} hours.', fontsize=16)
    ax[0, 0].set_title('Moisture Gas')
    ax[0, 1].set_title('Moisture Particle')
    ax[0, 2].set_title('Temp Gas')
    ax[0, 3].set_title('Temp Particle')

    temp_color = 'Reds'
    humidity_color = 'Blues'
    for index, h in enumerate(plot_times):
        ax[index, 0].set_ylabel(f'T: {int(h/(time-1) * hours)} hours                      ', rotation=0, size='large')
        im = ax[index, 0].imshow(moisture_gas_vector[h, :, :], vmin=moisture_gas_initial_bed,
                                 vmax=moisture_gas_initial_in, cmap=humidity_color)
        fig.colorbar(im, ax=ax[index, 0])

        ax[index, 1].imshow(moisture_particle_vector[h, :, :])
        im = ax[index, 1].imshow(moisture_particle_vector[h, :, :], vmin=moisture_particle_initial,
                                 vmax=moisture_particle_saturated, cmap=humidity_color)
        fig.colorbar(im, ax=ax[index, 1])


        ax[index, 2].imshow(temp_gas_vector[h, :, :])
        im = ax[index, 2].imshow(temp_gas_vector[h, :, :], vmin=temp_min, vmax=max_temp_gas, cmap=temp_color)
        fig.colorbar(im, ax=ax[index, 2])


        ax[index, 3].imshow(temp_particle_vector[h, :, :])
        im = ax[index, 3].imshow(temp_particle_vector[h, :, :], vmin=temp_min, vmax=max_temp_particle, cmap=temp_color)
        fig.colorbar(im, ax=ax[index, 3])
    plt.savefig('heatmap_chosen_times.pdf')
    # plt.show()

def slide_heat_map(
    moisture_gas_vector, moisture_particle_vector, temp_gas_vector, temp_particle_vector, temp_min, max_temp_particle,
    max_temp_gas, moisture_particle_initial, moisture_particle_saturated, moisture_gas_initial_bed, moisture_gas_initial_in, hours):
    # current layer index start with the first layer
    idx = 0
    time, height, length = np.shape(moisture_particle_vector)

    # figure axis setup
    fig, ax = plt.subplots(2, 2, figsize=(20, 13))
    fig.subplots_adjust(bottom=0.15)

    ax[0, 0].set_title('Moisture Gas')
    ax[0, 1].set_title('Moisture Particle')
    ax[1, 0].set_title('Temp Gas')
    ax[1, 1].set_title('Temp Particle')

    # display initial image
    im_y = ax[0, 0].imshow(moisture_gas_vector[idx, :, :], cmap='Blues', vmin=moisture_gas_initial_bed, vmax=moisture_gas_initial_in)
    fig.colorbar(im_y, ax=ax[0, 0])

    im_x = ax[0, 1].imshow(moisture_particle_vector[idx, :, :], cmap='Blues', vmin=moisture_particle_initial, vmax=moisture_particle_saturated)
    fig.colorbar(im_x, ax=ax[0, 1])

    im_tg = ax[1, 0].imshow(temp_gas_vector[idx, :, :], cmap='Reds', vmin=temp_min, vmax=max_temp_gas)
    fig.colorbar(im_tg, ax=ax[1, 0])

    im_tp = ax[1, 1].imshow(temp_particle_vector[idx, :, :], cmap='Reds', vmin=temp_min, vmax=max_temp_particle)
    fig.colorbar(im_tp, ax=ax[1, 1])

    # setup a slider axis and the Slider
    ax_depth = plt.axes([0.23, 0.02, 0.56, 0.04])
    slider_depth = Slider(ax_depth, 'Time', 0, time-1, valinit=0)

    fig.suptitle(f'Moisture & temperature in cylinder at time: {int(idx)} hours.', fontsize=16)
    # update the figure with a change on the slider
    def update_depth(val):
        idx = int(round(slider_depth.val))
        im_x.set_data(moisture_particle_vector[idx, :, :])
        im_y.set_data(moisture_gas_vector[idx, :, :])
        im_tp.set_data(temp_particle_vector[idx, :, :])
        im_tg.set_data(temp_gas_vector[idx, :, :])
        time_current = idx/time * hours
        fig.suptitle(f'Moisture & temperature in cylinder at time: {int(time_current)} hours.',
                     fontsize=16)

    slider_depth.on_changed(update_depth)
    plt.savefig('heatmap_slider.pdf')
    plt.show()
