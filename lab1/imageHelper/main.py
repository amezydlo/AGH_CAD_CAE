import tkinter as tk
from PIL import Image, ImageTk


def on_click(event):
    # Calculate coordinates on a 0-100 scale
    scale_x = (event.x / image_width) * 100
    scale_y = 100 - (event.y / image_height) * 100  # Invert Y scale
    clicks.append((scale_x, scale_y))
    print(f"Clicked on image at: ({scale_x:.2f}, {scale_y:.2f})")


def update_crosshair(event):
    # Clear canvas
    canvas.delete("crosshair")

    # Draw crosshair lines
    canvas.create_line(event.x, 0, event.x, image_height + margin, fill='red', tags="crosshair")  # Vertical line
    canvas.create_line(0, event.y, image_width, event.y, fill='red', tags="crosshair")  # Horizontal line


def draw_scale():
    # Draw scale on X axis
    for i in range(0, 101, 10):
        x = (i / 100) * image_width
        canvas.create_line(x + margin, image_height + margin, x + margin, image_height + margin + 10, fill='black')
        canvas.create_text(x + margin, image_height + margin + 20, text=str(i), anchor=tk.N, font=('Arial', 8))

    # Draw scale on Y axis (inverted)
    for i in range(0, 101, 10):
        y = (i / 100) * image_height
        canvas.create_line(margin - 10, image_height - y + margin, margin, image_height - y + margin, fill='black')
        canvas.create_text(margin - 20, image_height - y + margin, text=str(i), anchor=tk.E, font=('Arial', 8))


def on_closing():
    # Sort clicks by X coordinate and round to 2 decimal places
    sorted_clicks = sorted(clicks, key=lambda point: point[0])
    rounded_clicks = [(round(x, 2), round(y, 2)) for x, y in sorted_clicks]

    # Prepare data for spline_comb function
    x_values = [point[0] for point in rounded_clicks]
    y_values = [point[1] for point in rounded_clicks]

    # Repeat first and last X
    if x_values:
        first_x = x_values[0]
        last_x = x_values[-1]
        spline_comb_xs = [first_x] + x_values + [last_x]
    else:
        spline_comb_xs = []

    # Display in the desired format
    print(f"splines_comp(1000, {spline_comb_xs}, {y_values})")

    root.destroy()


# Create main window
root = tk.Tk()
root.title("Click on the image")

# Margin parameter
margin = 0

# Load image
image_path = "../image.jpg"  # Provide your image path here
image = Image.open(image_path)
image_width, image_height = image.size
photo = ImageTk.PhotoImage(image)

# Create canvas and add image
canvas = tk.Canvas(root, width=image_width + 2 * margin,
                   height=image_height + margin + 40)  # Extra space for the scale
canvas.pack()
canvas.create_image(margin, margin, anchor=tk.NW, image=photo)

# Hide cursor
canvas.config(cursor="none")

# List to store clicks
clicks = []

# Bind click and mouse movement events
canvas.bind("<Button-1>", on_click)
canvas.bind("<Motion>", update_crosshair)

# Draw scale
draw_scale()

# Set window close function
root.protocol("WM_DELETE_WINDOW", on_closing)

# Start main loop
root.mainloop()
