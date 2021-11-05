# Protein Analysis Using a Least-squares Approach
This is a Shiny app whose purpose is to measure the secondary structure composition of proteins from CD data. This was a tool originally developed for my Biophysical Chemistry class in late 2020. Since I found that it outperforms existing services with the same purpose, I decided to turn it into a web so that future students can obtain more accurate estimations. The original assignment, where the mathematical basis of this method is laid (which is little more than a non-linear least squares fit, really) can be found in the *www* folder in this repository. This tool was originally named SIMON, but since I made certain improvements I thought it would be better to give it a new name, an interface and leave it online for the world to use. Any questions regarding PAULA can be sent to my academic email address, malmriv at correo.ugr.es.

I wanted to make a Frequently Asked Questions section so that certain aspects of this service could be explained easily, but since no one as asked me a single question (I am writing this before launching the service) I will call it Reasonable Questions.

## **Reasonable Questions.**
### Why is the wavelength limited to the 170-250 nm range?
There are chemical and computational reasons for this. The chemical reason is that the interaction of far-UV light with proteins reveals the details concerning secondary information, so we need to work with a reasonable range of wavelengths. Infra-red light would be of no use for this purpose. The computational reason is that the idea underlying this service is just a non-linear least squares fit. The function being fitted is a combination of four interpolated spectra of single-structure proteins. That is, four spectra were experimentally measured by people smarter than me (Brahms, 1979), and I interpolated them and made a function that combines them in whatever proportion necessary. The original data taken by Brahms et al. is limited to that range, 170-250 nm. Therefore, it wouldn't make sense to compare data outside of this range with the original spectra.

### Why does it sometimes fail?
There are, again, several reasons. The first reason is that PAULA makes some simplifying assumptions about proteins. Even though there are many different patterns that can be observed in the secondary structure of a protein, I decided that four of them were sufficiently representative: beta sheets, beta turns, alpha helices and random coils. If you are dealing with a protein that contains some other structure, this service will force itself to find its spectrum as a combination of the four mentioned spectra, giving a bad result. The other reason is that non-linear fits and many computational methods usually require what's called a *seed*, a reasonable estimation from where to start. I wrote the program so that each execution uses a random but reasonable seed, so it would be a good idea to try several times if an initial attempt does not solve the spectra well enough. 

### Why the specific file formatting?
I took the files generated by the spectrometer which we used in the lab and tried to make this app so that almost no interaction is required from the student. It could have been a comma separated value, or an Excel file, but that would mean that future students will have to spend some time doing tedious work.

### Why aren't the units relevant when it comes to the amount of polarization?
The first thing I noticed when I started my weeks-long battle against the fourth assignment of this class was that no one really knew which units we were using. If you check online, you will find that at least four different units are currently in use. All of them can be converted between them. I found this is not important if we want to make a fit of the spectrum. The important thing is the shape of the spectrum, not *the size of it* (converting units would only scalate the spectrum by a certain amount). I took care of that by normalising the end result, obtaining a percentage for each different secondary structure.

### What's the meaning of life, the universe, and everything?
[42](https://en.wikipedia.org/wiki/Phrases_from_The_Hitchhiker%27s_Guide_to_the_Galaxy#The_Answer_to_the_Ultimate_Question_of_Life,_the_Universe,_and_Everything_is_42).	
