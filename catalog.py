# catalog.py
import math
from dataclasses import dataclass
from geometry import Vec3 

@dataclass
class CatalogStar:
    ra: float
    dec: float
    star_id: int     
    flag: str        
    magnitude: float 

    @property
    def spatial(self) -> Vec3:
        x = math.cos(self.dec) * math.cos(self.ra)
        y = math.cos(self.dec) * math.sin(self.ra)
        z = math.sin(self.dec)
        return Vec3(x, y, z)

def load_catalog(filepath: str) -> list[CatalogStar]:
    stars = []
    
    with open(filepath, 'r') as file:
        for line in file:
            if not line.strip():
                continue
                
            parts = line.split('|')
            
            if len(parts) == 5:
                ra_degrees = float(parts[0].strip())
                dec_degrees = float(parts[1].strip())
                
                star = CatalogStar(
                    ra=math.radians(ra_degrees),
                    dec=math.radians(dec_degrees),
                    star_id=int(parts[2].strip()),
                    flag=parts[3].strip(),
                    magnitude=float(parts[4].strip())*100 # why? its what lost does :shrug: i guess for speed?
                )
                stars.append(star)
                
    return stars