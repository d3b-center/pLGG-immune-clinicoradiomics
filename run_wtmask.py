import numpy as np
import nibabel as nb


def create_tumor_mask(img_path):
    """
    This function takes the tumor segmentation file (with different labels values: 1, 2, 3 and 4)
    as input and returns a binary whole tumor mask
    Parameters:
    img_path: This is the path to the input tumor segmentation file

    Returns:
    wt_mask: Binary whole tumor segmentation mask

    Example usage:
    wt_mask = create_tumor_mask('/path/to/input/seg/file')
    nib.save(qt_mask,'SegFileName.nii.gz')
    """
    # Load the tumor segmentation file
    img = nib.load(img_path)
    img_data = img.get_fdata()
    img_data = np.float32(img_data)

    # Create the binary mask by setting all label values greater than 0 to 1
    img_data[img_data > 0] = 1
    img_data = np.rint(img_data).astype('int16')

    wt_mask = nib.Nifti1Image(img_data, img.affine, img.header)

    return wt_mask
