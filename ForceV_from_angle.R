ForceV_from_angle <- function(target_angle = 5*pi/12, k, d=temp$d){
  
  tibble(H = sqrt(d^2 * (1 + tan(target_angle)^2)),
         ForceT = k*(H-d), 
         ForceV = ForceT*sin(target_angle)) %>%
    pull(ForceV)
  
}