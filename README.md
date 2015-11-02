So far not extensively tested.
I recommend testing on a copy of your flat Ntuples before writing output to the same file.

Let me know if you run into any problems.

## MonoJetCorrections

Either edit `default.cfg` to have the names of your directories, trees, and branches, or
you can make your own configuration file.
The corrector can then be run by doing
```
./ApplyCorrections.py -c <Config File Name>
```
If you leave out the `-c` option, `default.cfg` will be used.
Note that this will apply the footprint and smearing on all the events in the tree.
The values are only useful though in the appropriate control region selection.
Also remember, footprint corrections are only useful in data and smearing is only used in MC.

That's it! (hopefully)
