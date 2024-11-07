import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import pandas as pd

def calculate_framedrag_redshift(v_c, earth_radial_v=0):
    """Calculate redshift due to frame dragging"""
    radial_factor = 1 + earth_radial_v
    return radial_factor * (np.sqrt(1 / (1 - v_c**2)) - 1)

def calculate_doppler_redshift(v_c, earth_radial_v=0):
    """Calculate relativistic Doppler redshift"""
    v_total = (v_c + earth_radial_v)/(1 + v_c * earth_radial_v)
    return np.sqrt((1 + v_total)/(1 - v_total)) - 1

def calculate_combined_redshift(v_c, direction=1, earth_radial_v=0):
    """
    Calculate combined redshift with proper frame dragging dominance
    v_c: velocity as fraction of c
    direction: 1 for away, -1 for towards
    """
    # First calculate frame dragging effect
    gamma_fd = 1 / np.sqrt(1 - v_c**2)
    z_fd = gamma_fd - 1

    # Calculate effective velocity in the dragged frame
    v_eff = direction * v_c / gamma_fd

    # Calculate Doppler shift in dragged frame
    z_doppler = np.sqrt((1 + v_eff)/(1 - v_eff)) - 1

    # Combine effects (multiply wavelength ratios)
    total_ratio = (1 + z_fd) * (1 + z_doppler)
    z_total = total_ratio - 1

    return z_total

def generate_data_table(velocity, radius_parsecs, earth_v):
    """Generate data table for specific velocity value"""
    radius_factor = radius_parsecs / 100
    v_scaled = velocity * radius_factor

    data = {
        'Parameter': [
            'Angular Velocity (c)',
            'Scaled Velocity (c)',
            'Earth Radial Velocity (c)',
            'Distance (parsecs)',
            'Frame Dragging Redshift',
            'Doppler Away',
            'Doppler Towards',
            'Combined Away',
            'Combined Towards'
        ],
        'Value': [
            velocity,
            v_scaled,
            earth_v,
            radius_parsecs,
            calculate_framedrag_redshift(v_scaled, earth_v),
            calculate_doppler_redshift(v_scaled, earth_v),
            calculate_doppler_redshift(-v_scaled, earth_v),
            calculate_combined_redshift(v_scaled, 1, earth_v),
            calculate_combined_redshift(v_scaled, -1, earth_v)
        ]
    }
    return pd.DataFrame(data)

def update(val):
    """Update plot and data table"""
    radius_parsecs = radius_slider.val
    earth_v = earth_slider.val
    reference_parsecs = 100
    radius_factor = radius_parsecs / reference_parsecs

    # Scale velocities based on radius
    scaled_velocities = velocities * radius_factor

    # Mask out invalid velocities
    valid_mask = scaled_velocities < 0.99

    # Calculate redshifts
    framedrag = calculate_framedrag_redshift(scaled_velocities[valid_mask], earth_v)
    doppler_away = calculate_doppler_redshift(scaled_velocities[valid_mask], earth_v)
    doppler_towards = calculate_doppler_redshift(-scaled_velocities[valid_mask], earth_v)
    combined_away = calculate_combined_redshift(scaled_velocities[valid_mask], 1, earth_v)
    combined_towards = calculate_combined_redshift(scaled_velocities[valid_mask], -1, earth_v)

    # Update lines
    lines[0].set_data(velocities[valid_mask], framedrag)
    lines[1].set_data(velocities[valid_mask], doppler_away)
    lines[2].set_data(velocities[valid_mask], doppler_towards)
    lines[3].set_data(velocities[valid_mask], combined_away)
    lines[4].set_data(velocities[valid_mask], combined_towards)

    # Generate and print data table for 0.5c
    table = generate_data_table(0.5, radius_parsecs, earth_v)
    print("\nData at 0.5c reference point:")
    print(table.to_string(index=False))

    # Update title
    direction = "outward" if earth_v > 0 else "inward"
    ax.set_title(f'Redshift at {radius_parsecs:.1f} parsecs from axis\nEarth moving {direction} at {abs(earth_v):.2f}c', 
                fontsize=14)

    # Set fixed y-axis limits to avoid NaN issues
    ax.set_ylim(-2, 5)

    fig.canvas.draw_idle()

# Create figure and axis
fig, ax = plt.subplots(figsize=(12, 8))
plt.subplots_adjust(bottom=0.25)

# Create velocity array
velocities = np.linspace(0.01, 0.9, 1000)

# Initial plots with masked invalid velocities
valid_mask = velocities < 0.99
lines = [
    ax.plot(velocities[valid_mask], calculate_framedrag_redshift(velocities[valid_mask]), 'b-', 
           label='Frame Dragging Only', linewidth=2)[0],
    ax.plot(velocities[valid_mask], calculate_doppler_redshift(velocities[valid_mask]), 'r--',
           label='Doppler Only (away)', linewidth=1)[0],
    ax.plot(velocities[valid_mask], calculate_doppler_redshift(-velocities[valid_mask]), 'g--',
           label='Doppler Only (towards)', linewidth=1)[0],
    ax.plot(velocities[valid_mask], calculate_combined_redshift(velocities[valid_mask], 1), 'r-',
           label='Combined (away)', linewidth=2)[0],
    ax.plot(velocities[valid_mask], calculate_combined_redshift(velocities[valid_mask], -1), 'g-',
           label='Combined (towards)', linewidth=2)[0]
]

# Customize plot
ax.set_xlabel('Angular Velocity (fraction of c)', fontsize=12)
ax.set_ylabel('Redshift (z)', fontsize=12)
ax.grid(True, alpha=0.3)
ax.legend(fontsize=10)

# Add reference line
ax.axvline(x=0.5, color='gray', linestyle=':', alpha=0.5)
ax.text(0.51, ax.get_ylim()[0], '0.5c', fontsize=10, alpha=0.5)

# Create sliders
ax_radius = plt.axes([0.2, 0.1, 0.6, 0.03])
ax_earth = plt.axes([0.2, 0.05, 0.6, 0.03])

radius_slider = Slider(ax_radius, 'Distance from Axis (parsecs)', 10, 1000, valinit=100, valstep=1)
earth_slider = Slider(ax_earth, 'Earth Radial Velocity (c)', -0.5, 0.5, valinit=0, valstep=0.01)

# Connect sliders
radius_slider.on_changed(update)
earth_slider.on_changed(update)

# Initial update
update(None)

plt.show()