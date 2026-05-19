## get SST data and calculate DHW

library(tidyverse)
library(tidync) # For easily dealing with NetCDF data
library(rerddap) # For easily downloading subsets of data
library(doParallel) # For parallel processing
library(slider)
library(here)

# Bring in the in situ temperature data
Temp_LTER1<-read_csv(here("Data","raw_data","MCR_LTER02_BTM_Backreef_Forereef_20250522.csv")) %>%
  filter(reef_type_code == "Backreef") 

Mean_daily_Temp<-Temp_LTER1 %>%
  mutate(date = as.Date(time_local)) %>%
  group_by(date)%>%
  summarise(daily_temp =mean(temperature_c, na.rm = TRUE))


# The information for the NOAA OISST data
rerddap::info(datasetid = "ncdcOisst21Agg_LonPM180", url = "https://coastwatch.pfeg.noaa.gov/erddap/")


# Moorea NOrth shore 17°27'41"S 149°47'04"W
# -17.461389, 149.784444


# This function downloads and prepares data based on user provided start and end dates
OISST_sub_dl <- function(time_df){
  OISST_dat <- rerddap::griddap(datasetx = "ncdcOisst21Agg_LonPM180",
                                url = "https://coastwatch.pfeg.noaa.gov/erddap/", 
                                time = c(time_df$start, time_df$end), 
                                zlev = c(0, 0),
                                latitude = c(-17.46, -17.46),
                                longitude = c(149.78, 149.78),
                                fields = "sst")$data |> 
    dplyr::mutate(time = base::as.Date(stringr::str_remove(time, "T12:00:00Z"))) |> 
    dplyr::rename(t = time, temp = sst, lon = longitude, lat = latitude) |> 
    dplyr::select(lon, lat, t, temp) |> 
    stats::na.omit()
}

dl_years <- data.frame(date_index = 1:3,
                       start = c("2006-01-01", "2013-01-01", "2019-01-01" 
                                 ),
                       end = c("2012-12-31", "2018-12-31","2025-12-31" ))

# Download all of the data with one nested request
# The time this takes will vary greatly based on connection speed

  OISST_data <- dl_years |> 
    dplyr::group_by(date_index) |> 
    dplyr::group_modify(~OISST_sub_dl(.x)) |> 
    dplyr::ungroup() |> 
    dplyr::select(lon, lat, t, temp)
 # 1.19    0.19  360.18 

  OISST_data1<- OISST_data
   OISST_data1$lon <- as.vector(OISST_data$lon)
   OISST_data1$lat <- as.vector(OISST_data$lat)
  
  write_csv(x = OISST_data1, file = "Data/OISST_data.csv")

# plot it
#OISST_data |> 
#dplyr::filter(t == "2019-12-01") |> 
#  ggplot2::ggplot(aes(x = lon, y = lat)) +
#  ggplot2::geom_tile(aes(fill = temp)) +
#   ggplot2::borders() + # Activate this line to see the global map
#  ggplot2::scale_fill_viridis_c() +
#  ggplot2::coord_quickmap(expand = F) +
#  ggplot2::labs(x = NULL, y = NULL, fill = "SST (°C)") +
#  ggplot2::theme(legend.position = "bottom")


SST_year_month<-OISST_data %>%
  mutate(Year = year(t),
         Month = month(t)) %>%
  group_by(Year, Month)%>%
  summarise(mean_SST = mean(temp, na.rm = TRUE),
            max_SST = max(temp, na.rm = TRUE)) 

SST_year_month %>%
  mutate(Year_Month = paste(Year, Month))%>%
  ggplot(aes(x = Year_Month, y = max_SST))+
  geom_point()
  
SST_year<-OISST_data %>%
  mutate(Year = year(t),
         Month = month(t)) %>%
  group_by(Year)%>%
  summarise(mean_SST = mean(temp, na.rm = TRUE),
            max_SST = max(temp, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(Year_Benthic = Year + 1) # the benthic data should be last years temperature driving the change


SST_plot<-SST_year %>%
  filter(Year < 2026)%>%
  ggplot(aes(x = Year, y = mean_SST))+
  geom_point(size = 2)+
  geom_smooth(method = "lm", color = "black")+
  labs(x = "",
       y = " ")+
  theme_bw()+
  theme(axis.text = element_text(size = 20, face = "bold"),
        axis.title = element_text(size =18, face = "bold"),
        panel.grid = element_blank())
  
ggsave(here("Output","YearlySST.pdf"), width = 4, height = 4)
  

## Calcualte DHW by year
# ---- Helper to compute MMM from a chosen climatology window ----
compute_mmm <- function(df, date_col = "date", sst_col = "sst",
                        clim_years = 2006:2025) {
  df %>%
    transmute(date = .data[[date_col]],
              sst  = .data[[sst_col]],
              year = year(date),
              month = month(date)) %>%
    filter(year %in% clim_years) %>%
    group_by(year, month) %>%
    summarize(monthly_mean = mean(sst, na.rm = TRUE), .groups = "drop_last") %>%
    summarize(clim_monthly_mean = mean(monthly_mean, na.rm = TRUE), .groups = "drop") %>%
    summarize(MMM = max(clim_monthly_mean, na.rm = TRUE)) %>%
    pull(MMM)
}

# ---- Main function: compute daily DHW, then summarize by year ----
OISST_data<-OISST_data %>%
  rename(date = t, sst  = temp)

compute_dhw <- function(df,
                        date_col = "date",
                        sst_col  = "sst",
                        mmm = NULL,                 # set to numeric to override climatology
                        clim_years = 2006:2015,     # used if mmm is NULL
                        window_days = 84) {
  
  stopifnot(all(c(date_col, sst_col) %in% names(df)))
  dat <- df %>%
    transmute(date = as_date(.data[[date_col]]),
              sst  = as.numeric(.data[[sst_col]])) %>%
    arrange(date)
  
  # Optional: check daily regularity and warn if large gaps
  gaps <- diff(dat$date)
  if (length(gaps) && max(as.integer(gaps), na.rm = TRUE) > 3) {
    warning("Detected date gaps >3 days; DHW may be underestimated in gap periods.")
  }
  
  # MMM: compute if not provided
  if (is.null(mmm)) {
    mmm <- compute_mmm(dat, date_col = "date", sst_col = "sst", clim_years = clim_years)
  }
  
  # HotSpot (°C): full anomaly above MMM, but only days >= 1°C contribute to DHW
  # (NOAA CRW standard: accumulate the full hotspot value, not just the excess above 1°C)
  dat <- dat %>%
    mutate(
      hotspot     = pmax(0, sst - mmm),
      hotspot_crw = ifelse(hotspot >= 1, hotspot, 0)
    )
  
  # DHW daily (°C-weeks): trailing 84-day sum of qualifying hotspots / 7
  dat <- dat %>%
    mutate(
      dhw = slide_dbl(hotspot_crw,
                      .before = window_days - 1,
                      .complete = FALSE,
                      .f = ~ sum(.x, na.rm = TRUE) / 7)
    )
  
  # Annual summaries
  yearly <- dat %>%
    mutate(year = year(date)) %>%
    group_by(year) %>%
    summarize(
      dhw_max = max(dhw, na.rm = TRUE),                 # standard report metric
      dhw_cumulative = sum(dhw, na.rm = TRUE),          # total DHW across the year
      days_with_hotspot = sum(hotspot > 0, na.rm = TRUE),
      .groups = "drop"
    )
  
  list(
    mmm = mmm,
    daily = dat,     # columns: date, sst, hotspot, dhw
    yearly = yearly  # columns: year, dhw_max, dhw_cumulative, days_with_hotspot
  )
}

# Use satellite OISST data with the published NOAA CRW MMM for this pixel.
# MMM = 28.25°C is from the NOAA CRW 1985-1993 CoralTemp climatology (5km product)
# for the Moorea grid cell (-17.46°S, -149.78°W). We cannot derive this baseline
# from the available OISST data, which only starts in 2006 — too late and too warm.
# Source: https://coralreefwatch.noaa.gov/product/5km/
res <- compute_dhw(OISST_data, date_col = "date",
                   sst_col   = "sst",
                   mmm       = 28.25)
 res$mmm
 res$yearly
 head(res$daily)

res$daily %>%
  ggplot(aes(date, dhw)) +
  geom_line() +
  geom_area(fill = "firebrick", alpha = 0.6, color = NA) +
  labs(x = "", 
       y = "DHW (°C-weeks)", 
     #  title = "Mo'orea Daily Degree Heating Weeks (trailing 12 weeks)"
       )+
  #scale_fill_manual("firebrick")+
  theme_minimal()+
  theme(axis.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) 

ggsave(filename = here("Output","DHW2.png"), width = 8, height = 4)

write_csv(SST_year, here("Data","SST_year.csv"))

# Figure X. Daily Degree Heating Weeks (DHW) at Mo'orea, French Polynesia (17.46°S, 149.78°W), 2006–2025. DHW (°C-weeks) represents the accumulation of thermal stress over a trailing 12-week (84-day) window, calculated from NOAA Optimum Interpolation Sea Surface Temperature v2.1 (OISST; Huang et al. 2021; Huang et al. 2020) data. Following NOAA Coral Reef Watch (CRW) methodology, only days on which the sea surface temperature (SST) anomaly exceeded the Maximum Monthly Mean (MMM) by ≥1°C contributed to the running sum. The MMM climatology value of 28.25°C was derived from the NOAA CRW CoralTemp 1985–1993 baseline climatology for this grid cell. DHW values ≥4°C-weeks indicate likely coral bleaching; values ≥8°C-weeks indicate likely severe bleaching and mortality.

#References to include:
  
 # Huang, B., C. Liu, V. Banzon, E. Freeman, G. Graham, B. Hankins, T. Smith, and H.-M. Zhang, 2021: Improvements of the Daily Optimum Interpolation Sea Surface Temperature (DOISST) Version 2.1. Journal of Climate, 34, 2923–2939. https://doi.org/10.1175/JCLI-D-20-0166.1
#Huang, B., C. Liu, V.F. Banzon, E. Freeman, G. Graham, W. Hankins, T.M. Smith, and H.-M. Zhang, 2020: NOAA 0.25-degree Daily Optimum Interpolation Sea Surface Temperature (OISST), Version 2.1. NOAA National Centers for Environmental Information. https://doi.org/10.25921/RE9P-PT57_

## ---- Cross-check: Download official CRW DHW from ERDDAP ----
# Dataset: NOAA_DHW (CoralTemp 5km daily product)
# Moorea pixel: lat -17.50, lon -149.75 (nearest 0.05° grid node)
# Note: CRW uses a fixed 1985-1993 MMM climatology, so DHW values will differ
# from res$yearly (which uses OISST + the published MMM of 28.25°C). The key
# comparison is the *shape* of the time series and peak-year agreement.

crw_dl <- function(start, end) {
  rerddap::griddap(
    datasetx  = "NOAA_DHW",
    url       = "https://coastwatch.pfeg.noaa.gov/erddap/",
    time      = c(start, end),
    latitude  = c(-17.50, -17.50),
    longitude = c(-149.75, -149.75),
    fields    = c("CRW_DHW", "CRW_SST", "CRW_HOTSPOT")
  )$data |>
    dplyr::mutate(date = as.Date(stringr::str_remove(time, "T12:00:00Z"))) |>
    dplyr::select(date,
                  lon        = longitude,
                  lat        = latitude,
                  dhw_crw    = CRW_DHW,
                  sst_crw    = CRW_SST,
                  hotspot_crw = CRW_HOTSPOT)
}

# Download in two chunks to avoid ERDDAP timeout
crw_chunk1 <- crw_dl("2006-01-01", "2013-12-31")
crw_chunk2 <- crw_dl("2014-01-01", "2020-12-31")
crw_chunk3 <- crw_dl("2021-01-01", "2025-12-31")
crw_data   <- dplyr::bind_rows(crw_chunk1, crw_chunk2, crw_chunk3)

# Annual max DHW from CRW official product
crw_yearly <- crw_data |>
  dplyr::mutate(year = lubridate::year(date)) |>
  dplyr::group_by(year) |>
  dplyr::summarise(dhw_crw_max = sum(dhw_crw, na.rm = TRUE), .groups = "drop")

# Join with computed values and compare
dhw_compare <- dplyr::left_join(res$yearly, crw_yearly, by = "year") |>
  dplyr::select(year, dhw_computed = dhw_max, dhw_crw_max)

print(dhw_compare)

# Visual comparison
dhw_compare |>
  tidyr::pivot_longer(c(dhw_computed, dhw_crw_max),
                      names_to  = "source",
                      values_to = "dhw_max") |>
  dplyr::mutate(source = dplyr::recode(source,
                  dhw_computed = "Computed (OISST + MMM 28.25°C)",
                  dhw_crw_max  = "Official CRW (CoralTemp)")) |>
  ggplot2::ggplot(ggplot2::aes(x = year, y = dhw_max, color = source)) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::labs(x = "Year", y = "Annual Max DHW (°C-weeks)", color = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom")

# Save the CRW daily data for future use
write_csv(crw_data, here("Data", "CRW_DHW_Moorea.csv"))

