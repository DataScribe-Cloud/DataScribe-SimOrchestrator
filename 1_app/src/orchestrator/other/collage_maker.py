import os
from PIL import Image
import random

def create_image_collage(input_directory, output_path, collage_size=(3600*2, 3600*8)):
    
    # Get a list of all image files in the directory
    image_files = [f for f in sorted(os.listdir(input_directory)) if f.lower().endswith(('.png', '.jpg', '.jpeg', '.gif', '.bmp'))]
    
    # Shuffle the list to randomize the order
    #random.shuffle(image_files)
    #sort(image_files)

    # Create a blank canvas for the collage
    collage = Image.new('RGB', collage_size)

    # Calculate the number of rows and columns based on the number of images
    num_images = len(image_files)
    num_cols = int(num_images**0.5)
    num_rows = (num_images + num_cols - 1) // num_cols
    
    print('num_images=',num_images)
    print('num_cols and num_rows:',num_cols,num_rows)

    # Calculate the size of each cell in the collage
    cell_width = collage_size[0] // num_cols
    cell_height = collage_size[1] // num_rows

    # Paste each image onto the collage
    for i, image_file in enumerate(image_files):
        image_path = os.path.join(input_directory, image_file)
        img = Image.open(image_path)

        # Resize the image to fit the cell
        img = img.resize((cell_width, cell_height), Image.LANCZOS)

        # Calculate the position to paste the image
        row = i // num_cols
        col = i % num_cols
        x = col * cell_width
        y = row * cell_height

        # Paste the image onto the collage
        collage.paste(img, (x, y))

    # Save the collage
    collage.save(output_path)
    collage.show()

# Example usage:
input_directory = '../data/mics'
print(input_directory)
output_path = '../data/collage_t5000.jpg'
create_image_collage(input_directory, output_path)

