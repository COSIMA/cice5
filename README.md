## Overview
This repository contains the access/auscom fork of cice5 used in the ACCESS-ESM1.6 and ACCESS-OM2. It was forked from https://github.com/CICE-Consortium/CICE-svn-trunk/ which in turn captured the trunk from the subversion (svn) repository of the Los Alamos Sea Ice Model, CICE, including release tags through version 5.1.2.

ACCESS-ESM1.6 related code came from https://code.metoffice.gov.uk/trac/cice/browser/main/branches/pkg/Config/vn5.1.2_GSI8.1_package_branch/cice?order=name&rev=334#source (MOSRS account required)

More recent versions are found in the [CICE](https://github.com/CICE-Consortium/CICE) and [Icepack](https://github.com/CICE-Consortium/Icepack) repositories, which are maintained by the CICE Consortium.  

If you expect to make any changes to the code, we recommend that you work in the CICE and Icepack repositories.  Changes made to code in this repository will not be accepted, other than critical bug fixes.

There is [PDF documentation](https://github.com/ACCESS-NRI/cice5/blob/master/doc/cicedoc.pdf) available for CICE 5.1.2, however some changes were made to this fork to support coupling with ACCESS-OM2 and ACCESS-ESM1.6, Parallel IO, ERA5 Forcing, BGC modelling and for other updates. Some of these changes are described in the [ACCESS-OM2 Technical Report](https://github.com/COSIMA/ACCESS-OM2-1-025-010deg-report).

## Useful links
* **Wiki**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki

   Information about the CICE model prior to version 6 including how to obtain the code

* **Version Index**: https://github.com/CICE-Consortium/CICE-svn-trunk/wiki/CICE-Versions-Index-(older)

   Numbered CICE releases prior to version 6. 

* **Resource Index**: https://github.com/CICE-Consortium/About-Us/wiki/Resource-Index

   List of resources for information about the Consortium and its repositories as well as model documentation, testing, and development.
