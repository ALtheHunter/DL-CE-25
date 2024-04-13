function [rescaledImage, mmin, mmax] = rescaleImage(original_image)

rescaledImage = (original_image - min(original_image(:) ) ) ./ (max(original_image(:) ) - min(original_image(:) ) );
mmin = min(original_image(:));
mmax = max(original_image(:));
end