source("collect

## Gorilla-specific acceleration
gor.male.aggression = c('GO:0001820','GO:0002118','GO:0002124','GO:0002121')
gor.sex.dimorphism = c('GO:00030238','GO:0046661','GO:0008209','GO:00030521')
gor.change.diet = c('GO:0007586','GO:00044245','GO:0042476','GO:0048565')
gor.brain.size = c('GO:0007420','GO:0030900','GO:0021536','GO:0021871',
  'GO:0021895','GO:0021872','GO:0001890')
gor.immunity.genes = c('GO:0006955','GO:0002376','GO:0045087','GO:0045321','GO:0002250',
  'GO:0002520','GO:0050776','GO:0050778')
# Deceleration
gor.sperm.competition = c('GO:0046546','GO:0007286','GO:0030317','GO:0048232')


## Human-specific acceleration
hum.brain.size = c('GO:0007420','GO:0030900','GO:0021537','GO:0021871','GO:0021895',
  'GO:0021872','GO:0001890')
hum.cognition = c('GO:0050890')
hum.morphology = c('GO:0060348','GO:0001501')
hum.immunity.genes = gor.immunity.genes

## Gorilla-human overlap
gor.hum.brain.size = c('GO:0007420','GO:0030900','GO:0021872')

