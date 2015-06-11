Lineage animation R script

This is a quick hack to create a movie, showing an embryo movie image,
and the expression lineage tree at the same time.

It’s short enough that I think it’s simplest to just copy the script,
edit it, and run it (rather than trying to add command-line parsing.)

You will need several libraries installed; if they’re not installed,
type at the R prompt:

install.packages(c(“animation”, “png”, “tiff”))

Then, edit the configuration section. Hopefully, the options at the
start should be self-explanatory.

There are other things, such as the image size, movie brightness,
relative sizes of the images, etc., which could be tweaked by altering
the actual code.

Finally, run the script, by starting R in the same directory as the
script, and typing:

source(“lineageAnimation.r”)

The script is not very fast.

This should write a file called animation.gif in the current
directory. This can be opened in Fiji (or any other video-editing program,
presumably), and saved in your favorite format. (Saving directly as .mp4
doesn’t seem to work.) Mac users will probably want to save as QuickTime
(.mov) format.



