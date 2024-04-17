import numpy as np
import nibabel as nb


def normalize(img_path):
    """
    This function normalizes the input image intensity values between 0 and 255
    Parameters:
    img_path: This is the path to the input image to be normalized

    Returns:
    normalized_image: Image with intensity values normalized between 0 and 255

    Example usage:
    normalized_image = normalize('/path/to/input/image')
    nib.save(normalized_image,'FileName.nii.gz')
    """
    # Load the image
    img = nib.load(img_path)
    img_data = img.get_fdata()
    img_data = np.float32(img_data)

    # Get the 99.99th and 0th percentile values
    max_val = np.percentile(img_data, 99.99)  # Find the 99.99th percentile
    min_val = np.percentile(img_data, 0)

    img_data[img_data > max_val] = max_val

    # Normalize the intensity values between 0 and 255
    img_data = ((img_data - min_val) * 255.0) / (max_val - min_val)
    img_data = np.rint(img_data).astype('int16')

    normalized_image = nib.Nifti1Image(img_data, img.affine, img.header)

    return normalized_image
