This zip-file contains files used in paper:

Allen, Treb and Arkolakis, Costas. "Trade and the Topography of the Spatial Economy."
Quarterly Journal of Economics. 2014. Volume 129 (3), p.1085-1140.

The files are:

cfs_coor_revised.dat - an ASCII file containing the latitude and longitude of different CFS areas.
cfs_eth_revised.dat - an ASCII file containing data on the demographics of each CFS area. Source: IPUMS
cfs_trade_revised.dat - an ASCII file containing data on the value of bilateral trade flows by mode. Source: CFS.

cfs_est.m - a Matlab program that (1) estimates the travel time between each CFS area over the transportation network using the Fast Marching Method (FMM) algorithm and (2) estimates the cost of traveling over different types of transportation networks using a discrete choice framework along with the observed mode specific bilateral trade flows.

Note #1: Folder names will need to be altered in the program by each user.

Note #2: The program relies on the "Accurate Fast Marching Method" Matlab toolkit available at:

http://www.mathworks.com/matlabcentral/fileexchange/24531-accurate-fast-marching

fn_mode_111513.m - a Matlab program called by cfs_est_112513.m that implements the discrete choice estimation.

fn_trebreg_112513.m - a Matlab program called by cfs_est_112513.m that implements an OLS regression with two sets of fixed effects.

transportationnetworks_*_100113.png - a set of "speed" image files for rail, road, and water transportation networks in the United States. Please see Data Appendix B of the paper for information on the source of these files.

