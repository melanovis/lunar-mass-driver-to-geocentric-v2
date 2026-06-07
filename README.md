# lunar-mass-driver-to-geocentric-v2



<img width="1920" height="1080" alt="Image" src="https://github.com/user-attachments/assets/53b61856-a958-4223-b21c-2041cd78b08a" />

This was a project aimed at illuminating what Δv and T/W requirements are imposed on commercial rockets travelling from the moon to the earth. Specifically these are rockets which are launched from a lunar massdriver (MD) or some equivalent lifting vehicle. This is the third in two previous attempts to understand this problem which failed to consider various orbital effects such as lunar libration, the eccentricty of the lunar orbit and precession of the target orbit.

These are rockets which travel across cislunar space, executing multiple burns to eventually arrive at a circular parking orbit. After remaining there for some time they burn to rendezvous with a recieving station where their cargo payloads are seperated and the rockets are melted down into useful mass. This is the unidirectional MDR architecture, which is suspected as economically plausible due to the substantially reduced cost and complexity required for these lunar-ISRU derived rockets (likely Al+O2 redox) which do not need to negotiate earth's highly disagreeable atmosphere. These are more comparable to flying shipping containers than the fiery launch vehicles of our current period.

I wanted to show a continuous access massdriver (CAM) was possible, this is a system which can send material to some LEO orbit (where the first large commercial markets will be) no matter the orbital conditions such that mass-to-LEO throughput can be completely maximized (thus beating limited throughput alternatives such as OTVs). A railway in freefall.
The issue with this however is that orbital conditions are constantly changing, so to get a good picture into how these MDRs need to be scaled and how they should behave, we need to optimize multiple transits across the most influential orbital conditions. In this case it was determined that these should be the lunar orbital true anomaly (ν) and the longitude of the ascending node (Ω) of the LEO target orbit. For computational reasons I did not look to sweep across the lunar nodal precession and instead held it at the angle which would demand the most inclination change work to efficiently transit down to LEO.

These MDR transits are multiburn and follow elliptical restricted 3 body problem (ER3BP) dynamics both for the sake of simulation accuracy such that careful manoeuvring is available such that MDRs can preform expensive inclination changes at Δv-efficient altitudes and exploit effects like lunar gravity assists.

Here its suspected that 1200km as the minimum altitude range where debris mitigation becomes substantially simpler and less dramatic, with the target orbit being picked at an altitude of 1211.2 km which is semi-lunar synchronous to the lunar nodal precession cycle. To decide the altitude for the parking orbit above the 1211.2km receiving station orbit, I sampled various orbits for the time and Δv required to hohmann and rendezvous, making sure to account for the small but notable inclination difference incurred by the differences in precession rates between the orbits. We then decided to just draw a line at a flat 100m/s to get a 1371.2 km altitude parking orbit where we’d only have to wait a maximum of 57.36 hours for the phasing to be right for rendezvous.

<img width="7702" height="2871" alt="Image" src="https://github.com/user-attachments/assets/d19f597d-3646-4630-a700-5f8ec782e29e" />

In this simulation work the inclination (relative to earth's equator) of the LEO recieving station orbit was 28.59° which periodically aligns with the lunar orbital plane when its at Ω=90°. This periodic alignment opens the door to very substantial Δv savings but only for a short period, though definitely enough to bring down the mean. Here its additionally assumed the receiving station can change its inclination to 18.3 degrees very gradually over a 9.3 year period (to follow the lunar nodal precession cycle which occurs on an 18.6 year period due to the influence of the sun) which is considered doable due to the Δv provided to the station from unburnt propellent margins from the received MDRs. 
Just as an additional note, the precession of the target orbit is not caused by the influence of the sun but instead the earth's equatorial oblateness.

To begin designing a transit constellation, we need to look at the first part of their journey, that being the ejection from the lunar MD. 
A MD is an enormous asset, an investment comparable to a major continental railway. Ideally we should not be building multiple of them, in the bootstrapping period, instead we should be clever with the design and placement of a singular asset such that it can service the market regardless of the time. 
Informing this decision can only partially be done by orbital mechanics alone. Many more factors are required to justify the placement of an MD. But with that being noted, lets detour for a moment to discuss what parameters we have available to us for just the MD part of a MDR transit simulation. 
These include:
- Ejection velocity (v), the speed of material leaving the end of the MD relative to the stationary lunar surface.
- Orientation, or the azimuth (θ) and elevation (ε) we can point an MD at.
- Location, the latitude (φ) and longitude (λ) on the moon where we place the mouth of the MD. 
To gauge how an MDR transit constellation will behave for a single MD we must hold orientation and location fixed as the MD cannot change its aim or move along the lunar surface given its scale. But then how will we know if there are spots on the lunar surface which are preferred by the orbital mechanics of cislunar transiting? 
In order to answer this question we decided to split this simulation project into two phases. The first phase, which we’ll call the ‘allocation phase’ would have far less samples and leave these parameters free. We’ll then look at where the optimisation system places the MDs for each true anomaly-omega sample and if any patterns or rough regions emerge we’ll use those, along with further terrain feature justifications to pick a location and orientation for the more major ‘fixed phase’ where only the ejection velocity will be a free parameter.

Additional important constraints to note:
- Ejection velocity was constrained to a range of 2.32 - 3.5km/s.
- Longitude was constrained between -90 and 90 degrees, ie the earth-facing side of the moon so communications can reliably remain unbroken
- Azimuth was constrained between 0 and 180 degrees, meaning the MDs only fire retrograde in this phase
- Elevation was constrained to a maximum of 0.5 degrees
- No aerobraking is incorporated and MDRs will not reach below an altitude of 200km. This has been done as aerobraking in any context other than a surface-descent (and in this case it's a skim across the highest, thinnest part of the atmosphere) is almost impossible to accurately model. Furthermore there are political implications regarding what occurs if the MDRs accidentally deorbit which are unattractive.  This 200km collision boundary also applies to the moon, for which its unpreferable to get too close to due to the orbital perturbations potentially introduced by lunar masscons.

For the allocation phase, 64 samples (transits) were optimised, which revealled many bugs and had to be restarted twice. These 64 samples are for a range of 8 lunar ν values across 8 target orbit Ω values.
The resulting transits are shown all together in this animation, note the red vectors indicate the direction of burns.

<p align="center">
  <video src="https://github.com/user-attachments/assets/c887a5b7-3af5-41e3-af79-ee0dd0c395b2" width="100%" controls></video>
</p>

allocation phase results
<img width="5760" height="2975" alt="Image" src="https://github.com/user-attachments/assets/d9c5ac32-0736-4b53-a75b-013e4f76a621" />
<img width="5760" height="2975" alt="Image" src="https://github.com/user-attachments/assets/906670dc-3374-4a6a-b804-2e44176b55c7" />



prime abundance spots
<img width="4836" height="2068" alt="Image" src="https://github.com/user-attachments/assets/6cc8b7fe-d4f6-43fc-a582-8c58fdd56290" />


parameters
<img width="5760" height="2975" alt="Image" src="https://github.com/user-attachments/assets/f475c505-61a3-412f-962a-739d34018082" />

histogram
<img width="5760" height="2975" alt="Image" src="https://github.com/user-attachments/assets/c31bd526-2cba-4850-af04-ba3903861afd" />




No AI was used in the making of this project
