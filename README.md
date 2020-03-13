# How to populate the Milky Way with 2WDs+ from a Tres generated file

DWD_pop_BP_ugriz.C is the code produsing the distribution of DWDs in the Milky Way from a SeBa file.


## References

* spatial distribution of binaries: https://arxiv.org/pdf/astro-ph/0312193.pdf and https://arxiv.org/pdf/1806.03306.pdf
* star formation history: https://arxiv.org/pdf/astro-ph/9902148.pdf


To run the code:
```
./get_outputs
```

There are several additional options to give to the code:
* -g 25 sets the g magnitude limit to 25 (standard value = 23)
* -P 10 sets the maximum period to 10 days (standard value 100d)
* -R 0  sets the resolution for the output



