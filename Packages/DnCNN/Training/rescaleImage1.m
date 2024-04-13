function [rescaledImage, mmin, mmax] = rescaleImage1(original_image)

rescaledImage = original_image ./ max(original_image(:) );
mmin = min(original_image(:));
mmax = max(original_image(:));
end