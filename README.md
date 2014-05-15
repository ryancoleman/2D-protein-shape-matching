2D-protein-shape-matching
=========================
![![doi](https://zenodo.org/badge/3853/ryancoleman/2D-protein-shape-matching.png)](http://dx.doi.org/10.5281/zenodo.9874)

Did you ever want to find a protein that looks like a cat? Or a dinosaur? Well, now you can.

First, you'll need to make a bunch of precomputed images of proteins from many views. If you have a bunch of pdb files (one example in pdbs/966c.pdb), you can just run build_pdb_shape_database.py and it will generate them. It takes awhile, so I just put all the ones I generated in database/ I made them 256x256 one-bit images. Here is one:

![Random Protein](https://github.com/ryancoleman/2D-protein-shape-matching/blob/master/database/1a8i/1a8i.0.0.0.0.1.0.0.0.png?raw=true)

Now, you want to query them with something. The easiest way is to make a file like this:

```
001111110000000000
001111111000000000
001111111110000000
011111111110000000
111111111111100000
111111111111111110
011111111111111110
001111111111111100
001111111111111100
000011111111111000
000011111111111000
000011111111100000
000011111111100000
000111111111000000
000111111110000000
000111110000000000
001111000000000000
001110000000000000
001110000000000000
001100010000000000
000110000000000010
000000000000000000
```

put it in the shapes/ directory and run convert_image_shape.py to get this:

![South America](https://github.com/ryancoleman/2D-protein-shape-matching/blob/master/shapes/sa.png?raw=true)

Now you can run query_shape_database.py shapes/sa.png (or whatever your shape is called) and wait a long time while each image is compared using the Tanimoto coefficient. I didn't optimize the code because I don't plan on using this very much, it does take quite awhile (I built a super enormous database too), but then it shows you this:

![Best South American Protein](https://github.com/ryancoleman/2D-protein-shape-matching/blob/master/database/2ewb/2ewb.0.0.0.0.1.0.3.3.png?raw=true)

as the best match for South America in the database. Cool, right?

Depends on https://github.com/drj11/pypng, I just took the main routine and put it in png.py in the root directory, along with all the other code I wrote. Enjoy. Let me know if you do something cool. png.py is MIT, this code is all GPLv2. Go nuts.

Now you can also run query_shape_database.py shapes/sa.png top  to ignore the top half of the image. Other valid options are bottom, left or right.

Possible extensions
-------------------

Right now, the premade database is just rotations every 0.3 radians around each cardinal dimension. Could easily make a much bigger database with all possible combinations.

The code to do comparisons of the images is really slow as mentioned. Speeding it up is likely easy.

Obviously 2D projections of proteins aren't proteins, maybe somebody needs 3D shapes matched too. But lots of protein images are still made for two dimensional media.

