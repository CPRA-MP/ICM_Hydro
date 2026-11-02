Hydrology Subroutine
================
Stokka Brown, Brady Couvillion, Craig Conzelmann, Kim de Mutsert, Jordan
Fischbach, Christina Hunnicutt, Mark McKelvy, Philip Quibodeaux, Hugh
Roberts, Mallory Rodrigue, Jenni Schindler, Kevin Suir, Gordon Thomson,
Jenneke Visser, and Eric White

------------------------------------------------------------------------

Citation:

Brown, S., Couvillion, B., Conzelmann, C., de Mutsert, K., Fischbach,
J., Hunnicutt, C., McKelvy, M., Quibodeaux, P., Roberts, H., Rodrigue,
M., Schindler, J., Suir, K., Thomson, G., Visser, J., and White, E.
(2017). 2017 Coastal Master Plan: Appendix C: Modeling Chapter 3 -
Modeling Components and Overview. Version Final. (p. 2-7). Baton Rouge,
Louisiana: Coastal Protection and Restoration Authority.

------------------------------------------------------------------------

## Contents

1.  [HYDROLOGY SUBROUTINE](#10-hydrology-subroutine)
    1.  [MODEL DESIGN](#11-model-design)
    2.  [MODEL DOMAIN](#12-model-domain)
    3.  [COMPARTMENTS](#13-compartments)
    4.  [LINK NETWORK](#14-link-network)
    5.  [COMPARTMENT AND LINK
        ATTRIBUTES](#15-compartment-and-link-attributes)
    6.  [GOVERNING EQUATIONS](#16-governing-equations)
2.  [REFERENCES](#20-references)

## 1.0 Hydrology Subroutine

The 2017 Coastal Master Plan hydrology subroutine of the ICM is a mass
balance compartment model used to predict water level, salinity,
sediment, and other water quality constituents for 50 years into the
future. It is integrated with feedbacks to the morphology and vegetation
subroutine. The 2017 hydrology subroutine was built upon the 2012
Coastal Master Plan ecohydrology model. The spatial resolution of the
compartments was increased and new processes were developed and
included. The following is a brief summary of key changes since 2012.
For additional details regarding this subroutine and how it is
integrated into the ICM, refer to [Attachment C3-1: Sediment
Distribution](http://coastal.la.gov/wp-content/uploads/2017/04/Attachment-C3-1_FINAL_02.22.2017.pdf)
and [Attachment C3-22: ICM
Integration](http://coastal.la.gov/wp-content/uploads/2017/04/Attachment-C3-22_FINAL_03.07.2017.pdf).

### 1.1 Model Design

The 2012 eco-hydrology modeling used a mass balance approach whereby the
system is characterized using a set of compartments connected by links
that control the flow of constituents among the compartments.
Specifically, three mass balance, link node, compartment models (Meselhe
et al., 2013) were used covering different regions of the coast: (1)
Atchafalaya-Terrebonne (AA), (2) Chenier Plain (CP), and (3)
Pontchartrain-Barataria (PB; Figure 1). The PB model was developed in
the Formula Translating (Fortran) system programming language using the
multi-type compartment design that subdivides a hydrological compartment
into upland, marsh, and open water subcompartments. The AA and CP models
were developed in the Berkeley Madonna equation solver governed by the
same hydrodynamic equations but using the single-type compartment
design, whereby a compartment is composed of only a channel, marsh, or
open water rather than being subdivided (Meselhe et al., 2012). Figure 2
shows a schematic of the multi-type and single-type compartment designs.
Considering the timeline for model development for the 2012 coastal
master plan, existing models were used where available (e.g., PB).
Although there were differences among the models used for different
parts of the coast, they were standardized to the extent possible, and
any differences were determined to be negligible in terms of producing
50 year model outputs. The 2017 hydrology subroutine adopted the PB
Fortran multi-type compartment design. Fortran was used for the 2017
effort because it is a more flexible coding language with more
versatility for integrating model code when compared to Berkeley
Madonna. The PB multi-type compartment set up was chosen over that used
in the AA and CP regions of 2012 because PB was already coded in
Fortran.

The new hydrology subroutine was developed as a single coast wide model
(Figure 3) but can be subdivided into the AA, CP, and/or PB regions if
only one or two of the regions are required for a particular simulation.
Since the 2017 hydrology subroutine has adopted the PB region’s design,
only the differences in 2017 hydrology subroutine, compared to 2012
eco-hydrology model, are discussed in this overview.

![Satellite map of coastal Louisiana showing the 2012 Coastal Master
Plan eco-hydrology model compartments and domains, outlined in yellow
(Chenier Plain), red (Atchafalaya-Terrebonne), and blue
(Pontchartrain-Barataria)](https://github.com/CPRA-MP/ICM_Hydro/blob/ICM_Hydro_scg/ICM-Hydro_MP17_figs/ICM-Hydro_Fig3.png)  

<caption>

Figure 1. 2012 Coastal Master Plan eco-hydrology model compartments and
domains.
</caption>

![Diagrams of multi-type (Pontchartrain-Barataria) and single-type
(Atchafalaya-Terrebonne and Chenier Plain) compartment designs. Example
multi-type compartments are on the left and single-type compartments on
the right. Views include an aerial view and a cross-section view for
each compartment type. Multi-type compartments have more than one cell
type represented (e.g., upland, marsh, open water, or channel). Arrows
show flows into and out of the compartments such as
precipitation-evaporation, precipitation-evapotranspiration,
inter-compartment exchange, intra-cell exchange, tributary input, and
Mississippi
River.](https://github.com/CPRA-MP/ICM_Hydro/blob/ICM_Hydro_scg/ICM-Hydro_MP17_figs/ICM-Hydro_Fig1.png)  

<caption>

Figure 2. Multi-type (left) and single-type (right) compartment
</caption>

![Satellite map of coastal Louisiana showing the 2017 ICM Hydrology
Subroutine domain and compartments, outlined in
blue](https://github.com/CPRA-MP/ICM_Hydro/blob/ICM_Hydro_scg/ICM-Hydro_MP17_figs/ICM-Hydro_Fig2.png)  

<caption>

Figure 3. 2017 Coastal Master Plan Integrated Compartment Model (ICM) -
hydrology subroutine compartments and domain.
</caption>

### 1.2 Model Domain

For 2017, the southern boundary of the model domain was shifted further
offshore to approximately the 30 m contour to reduce variation in
salinity so a constant offshore boundary salinity could be applied and
to have all the regions the same “distance” offshore. The east and west
boundaries had negligible change, and the northern boundary in the AA
region was shifted north to include the entire Atchafalaya Basin area
between the East Atchafalaya Basin Levee and the West Mississippi River
and Tributaries Levee. See Figure 1 and Figure 3 for the spatial extent
of each model domain.

### 1.3 Compartments

Spatial resolution of the compartments in all three regions was
increased compared to the 2012 effort to better represent coastal
processes, such as sediment distribution (Table 1). This was made
possible due to advances in computational processing power. Uncertainty
in sediment distribution dynamics increases as the size of the hydrology
compartments increase. Excluding large offshore hydrology compartments
(which are included in Table 1), compartments specifically in the area
of the modeling domain in which most diversions are modeled were reduced
from an average size of 180 km<sup>2</sup> in the 2012 effort to 20
km<sup>2</sup> for the 2017 effort. It is in these compartments that
accurate sediment distribution is most critical for predicting project
benefits. The increased spatial resolution in these areas enhances
accuracy associated with projected diversion land area benefits.

Compartments in the AA and CP regions were delineated to convert them
from the single-type to the multi-type compartment design. Unlike the
2012 eco-hydrology model, the 2017 hydrology subroutine now has coast
wide coverage of multi-type compartments. Overall, there are a total of
946 compartments across the coast in the 2017 hydrology subroutine,
compared to 403 in the 2012 eco-hydrology model.

<caption>

Table 1. Summary of hydrology compartments per region and model; values
include all compartments in the model domain, including large offshore
compartments.
</caption>

|  | 2012 Eco-hydrology Model |  | 2017 ICM Hydrology Subroutine |  |
|----|----|----|----|----|
| Region | Total | Area Ranges in km<sup>2</sup> (average) | Total | Area Ranges in km<sup>2</sup> (average) |
| AA | 165 | 0.04 – 3361 (121) | 388 | 1.02 – 1477 (69) |
| CP | 149 | 0.44 – 1844 (86) | 245 | 0.40 – 2833 (123)<sup>1</sup> |
| PB | 89 | 2.2 – 5844 (716) | 313 | 1.84 – 3187 (160) |

<sup>1</sup> average area increased as a result of the additional large
offshore compartments for the CP region.

### 1.4 Link Network

The ICM’s hydrology subroutine utilizes a link-node type of hydraulic
exchange. Each node of the subroutine is a hydrologic compartment, which
consists of open water, marsh, and (optionally) upland areas. The links
represent the hydraulic pathways (bayous, bay inlets, canals, marsh
areas, etc.) connecting each compartment to neighboring compartments.
All modeled processes (rainfall/runoff, sediment distribution, water
quality sources/sinks, wave generation, etc.) take place within the
node/compartment, and the link network is used to convey flows and
constituents between compartments.

Both the 2012 eco-hydrology model and the 2017 hydrology subroutine link
networks consist of intracompartmental and intercompartmental link
types. The intracompartmental link types used in the 2017 hydrology
subroutine consist of rainfall/runoff processes in upland, marsh, and
open water areas, as well as exchange flow between the open water and
marsh areas within a compartment. Of the intercompartmental link types,
both the 2012 and 2017 versions include channel, lock, saltwater barrier
or tide gate, and pump link types. The 2017 hydrology subroutine
additionally includes weir; orifice, culvert, or bridge; overland marsh;
ridge or levee; and dormant link types. Dormant links are inactive links
in the model that become activated only after certain criteria are met
during the model simulation. Intercompartmental links connect the open
water portion of one compartment to the open water portion of a
neighboring compartment, with the exception of the overland marsh links.
Overland marsh links allow for overland flow across marsh areas, and
connect the marsh area of one compartment to the marsh area of a
neighboring compartment.

Although both 2012 and 2017 versions include lock link types, there was
an update to the lock operations in 2017. The 2012 eco-hydrology model
lock operations were either operated based on the recorded schedule, or
not operated and remained open 100% of the time. The 2017 hydrology
subroutine can operate locks based on differential stage, hourly
schedule, downstream stage, downstream salinity, or downstream stage and
salinity.

In addition to the hydraulic control rules that can now be simulated via
the aforementioned operational regimes; weirs, orifices, culverts, and
tide gates were also added to better represent the highly engineered
hydrology of Coastal Louisiana within the ICM’s hydrology subroutine.

Overland marsh and ridge or levee link types were added to the 2017
hydrology subroutine to allow for the propagation of extreme water
levels during a major riverine flood and/or tropical system’s storm
surge.

Dormant links allow the building of a delta or breaching of a barrier
island. For deltas, they were added to the fan-like compartment layouts
in areas of actively growing deltas and areas of potential delta growth
(e.g., from sediment diversion projects). As sediment is deposited and
accumulates on a compartment bed, the open water link invert elevation
increases, resulting in a decrease of the link’s flow capacity (i.e.,
cross-sectional area). Once the capacity of the open water link is
reduced to be equal to the capacity of the regime channel (i.e., the
minimum cross-sectional area required to facilitate the specified
diversion flow rate and particle size), the original open water link is
deactivated and the dormant link – with dimensions of the regime channel
– is activated. Additionally, sediment deposition will no longer occur
in this compartment and is instead pushed downstream to the next
compartment. For barrier islands, a dormant link is added to each island
to be activated if an island breach occurs, forming an inlet and
allowing exchange of water from offshore to the bay.

### 1.5 Compartment and Link Attributes

The 2012 eco-hydrology compartment bed elevations and link invert
elevations were calculated using the Digital Elevation Model (DEM)
developed by United States Geological Survey (USGS) for the 2012 Coastal
Master Plan. The 2017 hydrology subroutine compartment bed elevations
and link invert elevations were calculated using the updated 2017
Coastal Master Plan DEM. For additional information about this DEM,
please refer to [Attachment C3-27: Landscape
Data](http://coastal.la.gov/wp-content/uploads/2017/04/Attachment-C3-27_FINAL_03.10.2017.pdf).

For additional information regarding compartment delineations and
hydrodynamic links, refer to [Attachment C3-22: ICM
Integration](http://coastal.la.gov/wp-content/uploads/2017/04/Attachment-C3-22_FINAL_03.07.2017.pdf).

### 1.6 Governing Equations

The hydrodynamic formulations used in the 2017 hydrology subroutine are
unchanged from 2012 eco-hydrology, as these were based on well-known and
widely accepted hydrodynamic principles. These formulations can be found
in the 2012 Coastal Master Plan report Appendix D-1 (Meselhe et al.,
2012).

An extensive effort to improve the simulation of sediment distribution
among and within hydrology compartments was undertaken for 2017. In
2012, a single sediment accumulation value was calculated for each
compartment, and sediment was distributed within a compartment based on
a sediment distribution probability surface, which was based upon the
weighting of factors such as distance from sediment source, frequency of
inundation, and distance from water bodies. Also, the 2012 eco-hydrology
model assumed that open water bottoms (i.e., the bed) have an endless
supply of sediment for resuspension, while the 2017 hydrology subroutine
assumes a maximum TSS concentration and a limit to the amount of
erodible material available for resuspension. In both the 2012 PB
eco-hydrology model and the 2017 hydrology subroutine, the
intracompartment exchange between the open water and marsh allowed
sediment transfer to the marsh surface. In the 2017 modeling effort, an
additional source of sediment (from marsh edge erosion) is included in
the morphology subroutine; refer to the morphology section below for
additional information.

Both the 2012 and 2017 versions of the model code distribute sediment
throughout the model domain to account for tropical storms. The 2012
eco-hydrology model did not simulate tropical storms explicitly, and
therefore, to schematize the process of sediment distribution from
tropical events to the landscape, 1000 g/m<sup>2</sup>/year of sediment
was assumed delivered to each hydrology compartment (Couvillion et al.,
2013) for the duration of the model run. The 2017 hydrology subroutine
explicitly simulates tropical storm events by applying elevated water
levels (i.e., storm surge) at the offshore boundary and the storm’s wind
field temporally and spatially along the storm’s path. Therefore, the
sediment from tropical storms is delivered from the offshore
compartments to the marsh as sediment is resuspended from the bed due to
higher wave energy, and the marsh is inundated due to higher water
levels.

The 2017 Coastal Master Plan water quality constituents are predicted as
part of the hydrology subroutine. The predictions are based on
advection-diffusion equations of chemical species. Source/sink terms are
included to account for mass transfers due to chemical kinetic
processes. Water quality parameters, including total suspended solids
(TSS), salinity (SAL), water temperature (TMP), nitrate + nitrite
nitrogen (NO_3\_), ammonium nitrogen (NH_4\_), total inorganic
phosphorus (TIP), dissolved organic phosphorus (DOP), dissolved organic
nitrogen (DON), bluegreen algae (ALG), detritus (DET) are simulated. The
ICM water quality formulations are unchanged from the 2012 Coastal
Master Plan AA and CP regional models. These formulations can be found
in the 2012 Coastal Master Plan report Appendix D-1 (Meselhe et al.,
2012).

## 2.0 References

Couvillion, B. (2017). 2017 Coastal Master Plan Modeling: Attachment
C3-27: Landscape Data. Version Final. (pp. 1-84). Baton Rouge,
Louisiana: Coastal Protection and Restoration Authority.

Couvillion, B.R., Steyer, G.D., Wang, H., Beck, H.J. and Rybczyk, J.M.
(2013). Forecasting the effects of coastal protection and restoration
projects on wetland morphology in coastal Louisiana under multiple
environmental uncertainty scenarios. In: N. Peyronnin and D. Reed
(Eds.), Louisiana’s 2012 Coastal Master Plan Technical Analysis. Journal
of Coastal Research, Special Issue, No. 67 (pp. 29–50).

McCorquodale, A., Couvillion, B., Dortch, Freeman, A., M., Meselhe, E.,
Reed, D., Roth, B., Shelden, J., Snedden, G., Wang, H., and White, E.
(2017). 2017 Coastal Master Plan: Attachment C3-1: Sediment
Distribution. Version Final. (pp. 1-56). Baton Rouge, Louisiana: Coastal
Protection and Restoration Authority.

Meselhe, E., McCorquodale, J.A., Shelden, J., Dortch, M., Brown, T.S.,
Elkan, P., … Wang, Z. (2013). Ecohydrology component of Louisiana’s 2012
Coastal Master Plan: Mass-balance compartment model. Journal of Coastal
Research, 67, 16–28. <https://doi.org/10.2112/SI_67_2.1>

Meselhe, E., McCorquodale, J. A., Sheldon, J., Dortch, M., Brown, S.,
Davis, M., Wang, Z., Elkan, P. and Schindler, J. (2012). Eco-hydrology
model technical report (Appendix D-1). In Louisiana’s Comprehensive
Master Plan for a Sustainable Coast. Baton Rouge, LA: Coastal Protection
and Restoration Authority.

White, E.D., Meselhe, E, McCorquodale, A, Couvillion, B, Dong, Z,
Duke-Sylvester, S.M., & Wang, Y. (2017). 2017 Coastal Master Plan:
Attachment C3-22: Integrated Compartment Model (ICM) Development.
Version Final. (pp. 1-49). Baton Rouge, Louisiana: Coastal Protection
and Restoration Authority.
