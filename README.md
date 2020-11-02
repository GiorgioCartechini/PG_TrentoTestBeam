In the source code there is there is implemented the Sensitive detector scoring, but it does not work.
I used SteppingAction and EventAction classes, instead.

Geometry:
1) Absorber
2) ScintStart
3) Target
4) Lyso
5) NaI
6) ScintVeto1 for LYSO
7) ScintVeto2 for NaI

you can change the absorber material and size

you can change the target material between Yttrium (default) and copper and the target radius.
See macro run89Y_70.mac as example of 70 MeV proton beam on yttrium target. 
