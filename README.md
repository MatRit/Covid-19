## Covid-19 - A simple statistical model for predicting ICU load in exponential phases of the disease

**Matthias Ritter, Derek V.M. Ott, Friedemann Paul, John-Dylan Haynes, Kerstin Ritter**

**Abstract:** 
One major bottleneck in the ongoing COVID-19 pandemic is the limited number of critical care beds. Due to the dynamic development of infections and the time lag between when patients are infected and when a proportion of them enters an intensive care unit (ICU), the need for future intensive care can easily be underestimated. To infer future ICU load from reported infections, we suggest a simple statistical model that (1) accounts for time lags and (2) allows for making predictions depending on different future growth of infections. We have evaluated our model for public data from three regions, namely Berlin (Germany), Lombardy (Italy), and Madrid (Spain). For each region separately and before extensive containment measures made an impact, we first estimate the model parameters (i.e., ICU rate, time lag between positive testing and entering ICU, and average stay in ICU). Whereas for Berlin, an ICU rate of 6%, a time lag of 6 days, and an average stay of 12 days in ICU provide the best fit of the data, for Lombardy and Madrid the ICU rates were higher (18% and 15%), the time lag shorter (0 and 3 days), and the average stay in ICU shorter as well (4 and 8 days). These regional differences reflect quite reasonably the differing seriousness of the situation. The region specific models are then used to predict future ICU load assuming either a continued exponential phase with varying growth rates (0–-15%) or linear growth. The model can be adjusted as estimates of parameters develop and can thus also help to predict a potential exceedance of ICU capacity. Although our predictions are based on small data sets, disregard non-stationary dynamics, and have a number of assumptions, especially an exponential or linear development of cases, our model is simple, robust, and adaptable and can be updated when further data become available. 

The manuscript is currently under review. A previous version of this manuscript, where the model has only been applied to data from Berlin, can be found here: https://arxiv.org/ftp/arxiv/papers/2004/2004.03384.pdf


## Code structure
The main file is Covid19.m The script can be adapted to data from different cities, regions or countries. The respective data has to be saved in an additional xlsx sheet. Ranges for the parameters ICU rate, average time in ICU and time lag between reported infection and ICU admission have to be defined by the user. The model first estimates those parameters for the available data (by reducing the root mean squared error) and then predicts the future ICU load for a fixed time horizon (e.g., 2 months; can also be determined by the user) assuming different exponential growth rates or linear growth with different slopes. Optionally, capacity limits can be given. 

## Data
All data we used are public and can be downloaded from the *Berlin Senate Department for Health, Nursing and Equal Opportunities* (https://www.berlin.de/sen/gpg/service/presse/2020/} for Berlin, the *Presidenza del Consiglio dei Ministri -- Dipartimento della Protezione Civile* (https://github.com/pcm-dpc/covid-19/blob/master/schede-riepilogative/regioni/dpc-covid19-ita-scheda-regioni-20200421.pdf} for Lombardy, and *Datadista* (https://github.com/datadista/datasets/tree/master/COVID\%2019) for Madrid. We additionally provide here the xlsx sheets for Berlin, Lombardy and Madrid with data until April 20, 2020. 
