#---
# Practical 2 - Using R in Earth System modelling part II

# Tom Matthews (Original Authors: Murray Hoggett and Tom Pugh)
# 19.09.2012
#---


# Hello! This is the second R practical on working with climate data. In the last
# practical, we introduced the programming constructs of a "variable" and a
# "function". We then worked with 1 and 2 dimensional data to make some
# scatter plots. In this practical, we will introduce new functions for opening
# and working with a type of file called a NetCDF file, which contains 2
# dimensional and higher dimensional data. But don't panic, we will still just
# be storing this new data in variables, and working with the data with
# functions, so you already know all the theory. 

# Anyway, let's get started! 


## 2 (and 3) dimensional climate data - NetCDF files. ---
# As mentioned in the last R practical, the first thing we always need to do is
# set our working directory to the folder where we have saved the data.
# Otherwise R doesn't know where to look for the data.

# Set the working directory to the folder where the data is saved.

setwd("/home/jovyan/Earth_Systems_Practical2")

# You probably haven't heard of NetCDF as a file type. NetCDF files are files
# that normally end with the file extension .nc or sometimes .grd. NetCDF files
# are used for various types of climate data and geophysical data. The structure
# of data inside a NetCDF file is like a stack of spreadsheets on top of each
# other, where the x and y variables are normally latitude and longitude and the
# different layers in the stack are different time slices. The file will contain
# a 2-D or 3-D matrix which has the variables of interest and one vector for the
# values of latitude, another for the values of longitude and a third for time.
# The time axis doesn't always exist, in which case you have a single layer and
# no time vector

# To work with NetCDF files in R, you need to install a couple of packages. R is
# such a popular language for science and data analysis because it has thousands
# of packages which can be downloaded (for free!) with a single line of code.
# These packages range from forecasting time series for weather analysis or
# stock prices, to working with geospatial data, to higher performance
# parallelisation routines, to adding pictures of cats to plots (yes, really!).
# If you can think of a scientific problem, there's probably already a package for it in R.

# In R, you install a package once using the "install.packages()" function. This
# normally only needs to be done once on a computer. Then the package must be
# loaded into R with the "library()" function. This needs to be done every time
# you restart R. Let's load the ncdf4 package.

# Install the "ncdf4" package from the internet. This normally only needs to be
# done a single time on a computer, and then it is installed forever.

install.packages("ncdf4")

# You should see a load of both black and red text flash by on the screen. The
# important bit is the last few lines, which should say thing have been loaded
# OK.  The package has now been installed, but it still needs to be loaded into
# R (and loading the package needs to be done at the start of every script with
# a package in it). Let's do this now.

# Load the "ncdf4" package into R. This needs to be done at the start of every R
# script that uses the package.

library("ncdf4")

# You shouldn't see any output from running this line of code, but we can now use 
# the functions in the ncdf4 package. 

# What we are going to walk you through doing in the next few code chunks is the
# standard workflow for using a NetCDF file in R. This workflow is as follows:

# 1. Open the NetCDF file connection in R with the "nc_open()" function
# 2. Find out what is actually inside the NetCDF file with the "print()" function
# 3. Extract the variables we want to use from the NetCDF file using the "ncvar_get()" function
# 4. Close the NetCDF file connection in R with the "nc_close()" function.

# Each of the steps above is quite easy, so take it easy! 

# First, we want to open the NetCDF file and look at the ranges of different
# variables. Luckily this is quite easy in R. First we use the ncdf4 library to
# open a connection to the file and store this connection in a variable, and
# then we use the print function to print a summary of the file. Lets do this
# now. Remember if you haven't run the install.packages() and the library()
# functions, then R will complain that it "can't find" the function - which
# means R doesn't know what these functions are because they haven't been loaded
# yet. This is why we used library("ncdf4") above, to load them all in.

# Open a connection to the NetCDF file and store this connection in a variable
# called ncfile. (don't worry about what we mean by a "connection" to the file,
# this will become clear throughout the examples.). The data are Gross Primary
# Productivity (GPP) derived by upscaling observations from the current global
# network of eddy-covariance towers

ncfile <- nc_open("GPP_MR_1deg_annual_1982-2011_mean.nc") 

# Print the header of the NetCDF file (i.e. print the NetCDF file's metadata)

print(ncfile)

# Take a moment to read the text. Even though the output may look long and confusing (some
# files give much more than this), this is the same sort of output you will get
# from every NetCDF file you open with R, so it will become easy after you've
# done it a couple of times.

# The first few lines of the output say:
# 3 variables (excluding dimension variables):
#        double natural[lon,lat]   
#            units: fraction
#            missing_value: -9999

# This tells us there are three variables in this NetCDF file called: time_bnds,
# gpp and std. Some files may include many different variables. You will often
# see "double" and "float". This refers to how many decimal places the computer
# is storing numbers as. The important bit is that it is just a normal decimal
# number.

# The next paragraph of the print() output tells us information about the variables in the NetCDF file:
#      4 dimensions: ...

# This tells us that the dimensions of the file are latitude (stored in the
# variable called "lat"), longitude (stored in the variable called "lon"). There is also
# a time variable. And finally a forth variable called (e.g. bnds or level) - don't worry about this if you
# see it, we typically won't need them in this module.

# After this, we have a paragraph which starts with,
#  14 global attributes:

# In a NetCDF file, the global attributes are normally useful things like who
# made the file, when it was made, a description, and other useful information
# which is often called meta-data.

# We now understand a bit about our NetCDF file. The next step is to extract
# variables in the NetCDF file into variables in R. To do this we use the
# "ncvar_get()" function. The function arguments "nc=" refers to the NetCdf
# file, while the "varid=" function argument refers to the variable name *in the
# NetCDF* file. This is why we needed to use the print() function in the last
# step; to find out the variable names.

# Extract the 'lat' variable in the netcdf file, and store it in a variable called 'lat'.

lat <- ncvar_get(nc=ncfile, varid='lat')

# Extract the 'lon' variable in the netcdf file, and store it in a variable called 'lon'.

lon <- ncvar_get(nc=ncfile,varid='lon')

# Extract the 'gpp' variable in the netcdf file, and store it in a variable called 'gpp'.

gpp <- ncvar_get(nc=ncfile, varid='gpp')

# When this file was created, the data were flipped around (so the south pole is in the north),
# and we need to flip it back. You can ignore this really but just make sure to run the line!

gpp <- gpp[,180:1]


# We now have 3 variables in R holding the information we extracted from the
# variables inside the NetCDF file. We want to know what is inside these
# variables. We could just print their whole contents to the screen by just
# running their name in the usual way. However we know that some of these
# variables might be very, very big, and this might cause the computer to crash.
# A safer way is to use the "length()" function or the "dim()" function to see
# how long the variable is. Lets do this:

# Use the length function to see how long each variable is

length(lat)

length(lon)

length(gpp)

# OK, so that final variable gpp has 64800 entries. Probably a good thing
# we didn't just print that to the screen. Many netcdf files have millions of
# lines, especially if they include time as a variable.

# We have now extracted all the data we need from the NetCDF file connection. We
# can now close the connection to the file. This frees up some computer memory
# for us and helps stop confusion later. Note that we didn't need to do this in
# the last practical for the "read.csv()" and "read.table()" functions, which
# close the data file automatically.

# Close the NetCDF file connection

nc_close(ncfile)

# Now that we have extracted all the data in our NetCDF file, we are ready to start making some maps of the data! 

## Making our first map ---

# The data inside a NetCDF file (and therefore the data in our variables in R
# which we have extracted from the NetCDF files) is a rectangular grid of
# numbers called a matrix.

# Remember we extracted three data variables from the NetCDF file with the
# ncvar_get() function. Those variables were lat, lon and gpp.

# Let's remind ourselves of the dimensions of our gpp (gross primary productivity) data:

dim(gpp)

# This output means we have 360 rows and 180 columns in our matrix.
# If we check back to the length of the variables from the print() function we
# can see this corresponds to latitude and longitude

# We can easily make a picture from a matrix using the image.plot() function in the
# fields package. This function has an automatic colour bar. It also allows us to 
# specify the labels of the x and y axes 

# install the fields package, which as the image.plot() function. Remember, you only
# need to install the package the first time, after that you can just load it using the
# library() function. We will install the maps package, which has some useful functions for
# plotting maps

 install.packages("fields")
 
 install.packages("maps")
 
# Load the fields and maps packages with the library() function

library("fields")
 
library("maps")

# Make a plot of our matrix using the image.plot() function, which adds a
# colourbar automatically. The first argument is the values to use for the x
# axis, the second argument is the values for the y axis, and the third argument
# is the actual data values.

image.plot(lon, lat, gpp)

# Cool! There is also a useful function called map() in R
# which allows us to add country outlines. This is often a good idea to add, as
# it shows there is no data at all for Antarctica, Greenland, and the Sahara.
# Note the add=TRUE argument inside the map() function which adds to the current plot. 

# Same as the previous plot code but...

image.plot(lon, lat, gpp)

# ... add a low resolution set of country outlines over the top for context. 

map(database = 'world', lwd=1.5, add = TRUE, col='black')

# Hopefully you will have seen the country lines added. Try
# changing some of the other function arguments to customize the plot.

### Changing colour palettes

# There are a number of standard color palettes in R. The rainbow colour palette
# is standard in image.plot(). However a number of recent scientific articles
# have pointed out that rainbow colour palettes can draw the eye to unimportant
# features as some color changes act like contours differently to different
# people. Therefore simpler colour palettes can often be better for displaying
# scientific data. Let's install the RColorBrewer package to give us access to a
# load of useful color palettes. The RColorBrewer packages gives us access to
# the brewer.pal() function. This allows us to give a number of breaks in the
# colour palette we want (normally up to 9 or 10, depending on the palette) as
# the first function argument, and the name of the colour palette as the second
# argument. So brewer.pal(10, "RdBu") gives us a color palette from red to blue
# with 10 breaks in it. This is cool, but we know low numbers for our gpp
# data are normally displayed with cold colours in environmental sciences, so we
# use the rev() function to reverse the order, making us a blue to red colour
# palette. See below:

# Install the RColorBrewer package, to give us access to lots of extra colour
# paletts. 

 install.packages("RColorBrewer")
 
# Load the RColorBrewer into R with the library() function. 
 
library("RColorBrewer")
 
# The same plot and map commands as before, but with the brewer.pal(9,'RdBu)
# function, which gives us a Red to Blue color ramp with 9 levels. This is then
# enclosed with the rev() function, which reverses the order of the red-blue
# color ramp to give us a blue-red color ramp, which is more normal for
# Environmental sciences.

image.plot(lon, lat, gpp, col = rev(brewer.pal(9, "RdBu")))

map(database = 'world', add = T, lwd=1.5)

# You can find more about what color palettes are found in the RColorBrewer in
# this link
# [https://earlglynn.github.io/RNotes/package/RColorBrewer/index.html](https://earlglynn.github.io/RNotes/package/RColorBrewer/index.html).


### Code it yourself!

#Have a go at making the same plot (i.e. using image.plot() and map() but 
#with a yellow-green-blue color palette: "YlGnBu"). Also see if you can remember 
#from last week how to add in a title (call it whatever you want)

# The same image.plot() and map() code as before, but with a yellow-green-blue color palette.





# Masking NetCDF data

# As a next task, we will  make a difference map to show the changes between two
# different time slices. We will then use our subsetting skills to slice our
# data into regions of interest. We then cover how to use a mask in a separate
# NetCDF file to mask out areas and regions of arbitrary shape, allowing us to
# subset geographic regions such as individual continents, or individual habitat
# areas. We also show how to find summary statistics for any of these subsetted
# areas. Let's get started!

## Making a difference map ---

# A common thing to want to do is to see a change over time in some variable, or
# a difference between two variables. For example, how global temperature has
# changed over time, or how productivity has changed in London over time, etc.
# If you want to look at the difference in just one place, a line graph is a
# good choice (often called a "time series plot", as change over time is time
# series data). But what if you want to look at how the whole world has changed
# over time, or how two estimates of the same quantity differ for the globe? One
# good way of doing this is to look at the difference between arrays of data.
# Remember the "difference" just means to subtract one variable from the other.
# When dealing with differences over time it is often good to subtract the later
# time step from the earlier. Areas which are positive then show an increase
# over time in whatever variable we are considering, while negative values
# correspond to decreases over time in whatever variable we are considering.

# Let's try now to look at the difference between gpp and gpp_modis - we will use
# a difference map by subtracting one slice from another.

#First, we need to load in the new gpp_modis file
#(MODIS_GPP_2000-2009_mean_1deg.nc). This is GPP data derived from observations from the
#MODIS satellite. Ideally in science we want to have multiple different ways of
#deriving the same variable. If they give consistent results then that improves
#our confidence in those results.

# Open a file connection to the NetCDF file, extract variables and then close

ncfile_modis <- nc_open("MODIS_GPP_2000-2009_mean_1deg.nc") 

print(ncfile_modis)

lat_modis <- ncvar_get(ncfile_modis, 'lat')

lon_modis <- ncvar_get(ncfile_modis,'lon') 

gpp_modis <- ncvar_get(ncfile_modis, 'GPP')

nc_close(ncfile_modis)

# Looking at the output of the print() function on our modis file reveals that the 
# GPP variable is in units: kg C m-2 yr-1. However, in our original gpp variable
# from the file before, you might remember the units were in: kg m-2 s-1. Note
# the differences in terms of yr-1 vs.s-1. To deal with this, lets correct our
# original gpp variable so that it is also in units of yr-1, which is conceptually easier
# to deal with when looking across annual timescales.

sec_in_day <- 86400 #Number of seconds per day

sec_in_yr <- sec_in_day * 365 #multiply by number of days in year

gpp_yr <- gpp * sec_in_yr

# Now we can try and make the plot! First, we need to find the difference
# between the two matrices of data, and we need to remember to use our new
# gpp_yr variable to ensure the units are the same in both variables

gpp_difference <- gpp_yr - gpp_modis

# Then we can make a plot of the difference

image.plot(gpp_difference)

#QUESTION: can you interpret what this plot is showing?

## Summary statistics on data ---

# You might well want to take summary statistics for a grid of data. This is
# incredibly easy in R. Say, for example, you wanted to find the mean of our
# gpp_difference, i.e. we wanted to find the average global difference between
# gpp_yr and gpp_modis:

# Find the mean of the differences between our two slices. We need to remember to use the na.rm = TRUE,
# argument as there are NAs in our gpp_difference file. These come from, in this case,
# missing data, i.e. regions of the world where data were not available

mean(gpp_difference, na.rm=T)

# Similarly, lets find some other summary statistics:

median(gpp_difference, na.rm=T)

sd(gpp_difference, na.rm=T) # the sd function returns the standard deviation. 

max(gpp_difference, na.rm=T)


## Subset only equatorial regions---

# You might want to find the mean (average) GPP over a smaller area than the whole world. This is easy to do. 

# Make a temporary variable with the original gpp data and 
# correct the units again at the same time.

temp_gpp <- gpp * sec_in_yr

# We subset our grid with the square bracket notation. The part within the square 
# brackets before the comma is used to select rows (longitude values here) and the 
# part after the comma is used to select columns (latitude values), e.g: region <- temp_gpp[LON, LAT]

# We can combine the use of longitude and latitude values to select a particular region.
# For example, lets subset both latitude and longitude to find Australia. 

australia_gpp <- temp_gpp[290:335 , 45:80]

# Make a quick plot of the data.

image.plot(australia_gpp)

# You might be thinking that those longitude and latitude values are not correct, i.e. they
# don't match with a true world map. Lets have another look at one

image.plot(lon, lat, gpp)

# And yeah Australia is between roughly 100 and 150 longitude rather than 290:335. This is
# because our temp_gpp matrix starts at 0 at the bottom corner and goes up
# to 180 at the top corner (rather than -90 to +90), and then from 0 at the left to 360 
# on the right (rather than -180 to + 180) - this is because you can't have a negative number as a row
# or column number (e.g. try and extract the -5th row of a table). So, when subsetting the matrix, we
# need to use 0 for -180 longitude and -90 latitude, and work it out from there

### Code it yourself!

# Using the same approach, try subsetting and plotting Europe. 
# Use image.plot(lon, lat, gpp) to help you work out the correct LON and LAT values for Europe,
# but remember that x-axis values you want go from 0 to 360, while the y-axis values go from
# 0 to 180.
# Try giving it a title, new axis titles, and changing the colour pallete



# We can also calculate summary statistics using just our subsetted region, e.g.

mean(australia_gpp, na.rm = TRUE) #remember to set na.rm = TRUE to remove squares with NA values.

### Code it yourself.

#calculate the mean for Europe and compare to Australia



## Subset using a mask ---

# So we can now subset any type of rectangular box or slice we want from out
# matrix of data that we extracted from a NetCDF file. But what if we want an
# area more complex than a rectangular box? This is best done with a mask. A
# mask is normally a NetCDF file filled with the number 1 everywhere in a region
# we want to extract, and filled with NA everywhere else. If we extract our mask
# into a variable, and then multiply our mask by our data, then everywhere
# multiplied by 1 will stay the same, while everywhere filled with NA will
# become NA. This is because in R, anything multiplied by NA becomes NA. (This
# is kinda like how in maths, anything multiplied by 0 becomes 0. But 0 can
# actually have meaning in climate data, so we don't want to set everywhere
# outside our mask to zero.).

# Let's quickly see an example. This example is sort of like marking out one
# place we want to get rid of:

# Demonstrating the idea of how to subset data in R using a mask. We want to
# have a mask filled with NAs where we want to not have data, and 1s where we do
# have data.

pretend_data1 <- c(0,1,2,3,4,5)

pretend_mask1 <- c(1,1,NA,1,1,1)

# When we multiple our mask and our data, the data will stay the same where it
# is being multiplied by 1, but will change to NA where it is multiplied by NA.
# This example demonstrates masking out just one value.

pretend_data1 * pretend_mask1

# Hopefully the idea of masking to get rid of data is pretty clear! Now lets
# actually do it with a real mask. Here we skip a load of steps in reading the
# NetCDF as I already know what is in the file. I still print the header of the
# NetCDF file, as there is almost always some useful information in the header:

# Read in a NetCDF file containing a mask for landcover types around the world.
# Here we use the same NetCDF workflow as shown above, but much reduced.
# Open a file connection to the NetCDF containing the mask values

ncfile = nc_open("ESACCI-LC-L4-LCCS-Map-300m-P1Y-2015-v2.0.7_1deg.nc")

# Print the header information for the NetCDF file

print(ncfile)

# OK, so immediately we have found useful information in the header. The header says:
# >1 variables (excluding dimension variables):
#         double landcover[lon,lat]   
#             0: No data
#             10: Cropland
#             11: Herbaceous cover
#             20: Cropland, irrigated or post-flooding
#             ... and so on...

# So it seems like this file is a mask for many different regions (in this case
# landcover types). Each region seems to be filled with a common value, going
# from 1 to 220. Let's try extracting the values into a variable called "mask"
# and plotting the mask to see if we are correct:

# Use the ncvar_get() function to get the values of the mask

mask <- ncvar_get(ncfile, 'landcover')

# Close the connection to the NetCDF file.

nc_close(ncfile)

# Make a quick plot of the values in the mask file 

image.plot(lon, lat, mask, col=rainbow(36), legend.lab="Number of the zone in the mask")

map(database = 'world', add = T, lwd=1.5)

# Yep, looks like we are correct. However, as we talked about earlier, we
# preferably want this data mask to be in the form of 1s where we want to
# extract the data, and NAs everywhere else. Luckily we can do this in 1 line of
# code in R. Meet the "ifelse()" function. The ifelse() function takes a
# condition as it's first argument (for example mask==3, which means select only
# the places where the mask equals 3), the value to change to when that
# condition is true as it's second argument (which will always be 1 for us), and
# the value to change to when that condition is false as it's third argument
# (which will always be NA for us). Let's see this in action to make a mask for
# cropland, select those by multiplying our data by our new mask, and then
# making a map of our selected data:

# Here we make a new mask for cropland using the ifelse() function. This line of
# code makes a matrix where cropland areas are set to 1, while the rest of the
# matrix is set to NA. This matrix is then saved to the variable cropland_mask.

cropland_mask <- ifelse(mask == 10, 1, NA)

# Select that data by multiplying the data in gpp_yr by the mask. 

gpp_yr_cropland <- gpp_yr * cropland_mask

# Make a plot of our subsetted data. 

image.plot(lon, lat, gpp_yr_cropland, legend.lab="GPP (kg C m-2 yr-1)")

map(database = 'world', add = T, lwd=1.5)

# That's all fine. But you'll notice that there are a lot of very similar
# landcover classes in that file, including several including cropland. e.g.
# 10: Cropland
# 11: Cropland, rainfed, herbaceous cover
# 12: Cropland, rainfed, tree or shrub cover
# 20: Cropland, irrigated or post-flooding
# So what if you wanted to, for example, make a mask for all cropland classes? In this
# case we can make selections for each class and add them all together. Here we
# need to use a common coding operator called "or", which is written "|". You can chain
# multiple "or" together within an ifelse() statement. Below we have chained
# four, one for each main cropland class. Logical statements like this are very
# powerful, together with "and" (&) you can test for an enormous number of
# conditions.

allcropland_mask <- ifelse(mask == 10 | mask == 11 | mask == 12 | mask == 20, 1, NA)

# Now use this modified mask to select data in gpp_yr

gpp_yr_allcropland = gpp_yr * allcropland_mask

# Make a plot of our subsetted data. 

image.plot(lon, lat, gpp_yr_allcropland, legend.lab="GPP (kg C m-2 yr-1)")

map(database = 'world', add = T, lwd=1.5)

#QUESTION: can you interpret this plot?

# And then after we have subsetted/masked/selected (these words all mean the
# same thing here) our data, we can find summary statistics for it. In case it
# has become confusing why we are finding summary statistics of things, it's
# because often a summary statistic is the best way of answering a scientific
# question. Obviously a scientific question isn't normally "can you subset
# cropland?" - which after doing the tutorial, the answer is yes, you can subset
# cropland! The scientific question is more likely to be something along the
# lines of "What is the average GPP of global croplands". To answer this, you
# would want to take the mean (another way of saying the average) of the
# croplands. So let's do this, lets find  Don't forget the na.rm=T (or
# na.rm=TRUE) function argument!

# Summary statistics are things like mean, median, maximum, minimum,
# interquartile ranges, etc. Summary statistics are often the best way to answer
# scientific questions, such as "what is the maximum modeled temperature in a
# tropic broadleaf forest". Let's see an example of how easy it is to get a
# summary statistic for the average yield of cropland areas

mean(gpp_yr_allcropland, na.rm = T)


###--- And we are finished! ---
#As with the first tutorial, I strongly recommend that you go back through
#this tutorial at least one more time. There is just too much information to
#take in in one read through. But seriously, nice work! You can now do a huge
#amount of things in R. You can manipulate data to change it's units via simple
#mathematical operations. You can install and load new packages for free. You
#can read and work with data files hundreds of times too big to open in excel.
#You can make and save publication quality maps and graphs. You can select, mask
#and filter data you want. You can find summary statistics for areas and regions
#of the world, and you can make several types of plot to summarize it. Well
#done! Good luck in using these skills for your own work, and enjoy!

#start trying to play with the code. Change bits on each line to see
#what happens. Get used to R throwing error messages, and get used to googling
#how to do things. From reading these guides you will now have a good idea of
#the sorts of things you can do with R.





