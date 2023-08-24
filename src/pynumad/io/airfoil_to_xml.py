def airfoil_to_xml(coords, reftext, fname):
    """WriteNuMADAirfoil  Write NuMAD airfoil files

    Parameters
    ----------
    coords : array
        Nx2 array of airfoil coordinate data.  First column contains
        x-values, second column contains y-values.  Airfoil coordinates are in
        order as specified by NuMAD (i.e. trailing edge = (1,0) and leading
        edge = (0,0)
    reftext : string
        string representing reference text
    fname : string
        full filename, incl extension, of NuMAD airfoil file to write

    Returns
    -------
    None
    """
    with open(fname, "wt") as fid:
        fid.write("<reference>\n%s</reference>\n" % (reftext))
        fid.write("<coords>\n" % ())
        for i in range(coords.shape[0]):
            fid.write("%8.12f\t%8.12f\n" % tuple(coords[i, :]))
        fid.write("</coords>" % ())
