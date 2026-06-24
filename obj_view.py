import vtk
import os
import glob
import numpy as np

def napravi_obj_animaciju_sa_koliderom(folder_path, fps=10):
    """
    Učitava sve .obj fajlove iz foldera i prikazuje ih kao animaciju.
    Dodaje crvenu sferu kao kolider na poziciji (3, 2, 3) sa prečnikom 0.7.
    
    Args:
        folder_path (str): Putanja do foldera sa .obj fajlovima
        fps (int): Broj frejmova u sekundi
    """
    
    # Pronađi sve .obj fajlove u folderu
    obj_fajlovi = sorted(glob.glob(os.path.join(folder_path, "*.obj")))
    
    if not obj_fajlovi:
        print("Nema .obj fajlova u folderu!")
        return
    
    print(f"Pronađeno {len(obj_fajlovi)} .obj fajlova")
    
    # Kreiraj renderer, render window i interactor
    renderer = vtk.vtkRenderer()
    render_window = vtk.vtkRenderWindow()
    render_window.AddRenderer(renderer)
    render_window.SetSize(800, 600)
    render_window.SetWindowName(".obj Animacija sa koliderom - VTK")
    
    interactor = vtk.vtkRenderWindowInteractor()
    interactor.SetRenderWindow(render_window)
    trackball_style = vtk.vtkInteractorStyleTrackballCamera()
    interactor.SetInteractorStyle(trackball_style)

    # Podesi pozadinu
    renderer.SetBackground(0.1, 0.2, 0.4)  # Tamno plava
    
    # Podesi osvetljenje
    light = vtk.vtkLight()
    light.SetPosition(1, 1, 1)
    renderer.AddLight(light)
    
    # Dodaj dodatno osvetljenje da se bolje vidi sfera
    light2 = vtk.vtkLight()
    light2.SetPosition(-1, 1, 1)
    renderer.AddLight(light2)
     
    # =========================================================
    # KREIRANJE domena (Box)
    # =========================================================

    box_source = vtk.vtkCubeSource()
    box_source.SetCenter(2.0, 2.0, 2.0)      # Pozicija centra (x, y, z)
    box_source.SetXLength(4.0)               # Dužina po X osi
    box_source.SetYLength(4.0)               # Dužina po Y osi
    box_source.SetZLength(4.0)               # Dužina po Z osi
    box_source.Update()

    # Kreiraj mapper za Box
    box_mapper = vtk.vtkPolyDataMapper()
    box_mapper.SetInputConnection(box_source.GetOutputPort())

    # Kreiraj actor za Box
    box_actor = vtk.vtkActor()
    box_actor.SetMapper(box_mapper)

    # Podesi boju - PLAVA (možete promeniti RGB vrednosti po želji)
    box_actor.GetProperty().SetColor(0.8, 0.8, 1.8)

    # 🔑 Ključna linija za providnost (0.0 = potpuno providno, 1.0 = neprovidno)
    box_actor.GetProperty().SetOpacity(0.05)

    # Opciono: Podesi osvetljenje za realističniji providni izgled
    box_actor.GetProperty().SetAmbient(0.2)
    box_actor.GetProperty().SetDiffuse(0.8)
    box_actor.GetProperty().SetSpecular(0.1)      # Manje sjaja da ne bi "blještalo" kroz providnost
    box_actor.GetProperty().SetSpecularPower(10)

    # Dodaj Box u renderer
    renderer.AddActor(box_actor)   
    
    # =========================================================
    # KREIRANJE CRVENE SFERE (KOLIDER)
    # =========================================================
    print("Dodavanje crvene sfere kao kolider...")
    
    # Kreiraj sferu
    sfera_source = vtk.vtkSphereSource()
    sfera_source.SetCenter(2.0, 1.2, 2.0)  # Pozicija (x, y, z)
    sfera_source.SetRadius(0.25)  # poluprečnik 0.25
    sfera_source.SetThetaResolution(32)  # Rezolucija sfere (što veće, to glatkija)
    sfera_source.SetPhiResolution(32)
    sfera_source.Update()
    
    # Kreiraj mapper za sferu
    sfera_mapper = vtk.vtkPolyDataMapper()
    sfera_mapper.SetInputConnection(sfera_source.GetOutputPort())
    
    # Kreiraj actor za sferu
    sfera_actor = vtk.vtkActor()
    sfera_actor.SetMapper(sfera_mapper)
    
    # Podesi boju - CRVENA (neprovidna)
    sfera_actor.GetProperty().SetColor(1.0, 0.0, 0.0)  # RGB: crvena
    sfera_actor.GetProperty().SetOpacity(1.0)  # Neprovidna
    
    # Dodaj malo sjaja da bude vidljivija
    sfera_actor.GetProperty().SetAmbient(0.3)
    sfera_actor.GetProperty().SetDiffuse(0.8)
    sfera_actor.GetProperty().SetSpecular(0.5)
    sfera_actor.GetProperty().SetSpecularPower(20)
    
    # Dodaj sferu u renderer
    renderer.AddActor(sfera_actor)
    
    # =========================================================
    
    # Lista za čuvanje svih aktora (modela fluida)
    fluid_aktori = []
    
    # Učitaj sve .obj fajlove (fluid)
    for i, obj_fajl in enumerate(obj_fajlovi):
        print(f"Učitavanje {i+1}/{len(obj_fajlovi)}: {os.path.basename(obj_fajl)}")
        
        # Čitač za .obj fajl
        reader = vtk.vtkOBJReader()
        reader.SetFileName(obj_fajl)
        reader.Update()
        
        # Proveri da li je učitavanje uspelo
        if reader.GetOutput().GetNumberOfPoints() == 0:
            print(f"Upozorenje: {os.path.basename(obj_fajl)} nema tačaka ili nije ispravan .obj fajl")
            continue
        
        # Kreiraj mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())
        
        # Kreiraj actor
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        
        # Podesi boju za fluid (plava, providna da se vidi sfera kroz nju)
        actor.GetProperty().SetColor(0.0, 0.5, 1.0)  # Plava
        actor.GetProperty().SetOpacity(0.7)  # Delimično providna (0 = potpuno providno, 1 = potpuno neprovidno)
        actor.GetProperty().SetAmbient(0.2)
        actor.GetProperty().SetDiffuse(0.8)
        actor.GetProperty().SetSpecular(0.3)
        
        # Sakrij sve aktore fluida osim prvog
        if i > 0:
            actor.VisibilityOff()
        
        renderer.AddActor(actor)
        fluid_aktori.append(actor)
    
    if not fluid_aktori:
        print("Nijedan .obj fajl nije uspešno učitan!")
        return
    
    print(f"Uspešno učitan {len(fluid_aktori)} modela fluida")
    print(f"Crvena sfera (kolider) dodata na poziciji (3, 2, 3) sa prečnikom 0.7")
    
    # Podesi kameru da vidi celu scenu
    renderer.ResetCamera()
    
    # Pomeri kameru malo da bolje vidi sferu
    kamera = renderer.GetActiveCamera()
    kamera.SetPosition(10, 5, 10)  # Pozicija kamere
    kamera.SetFocalPoint(3, 2, 3)  # Gledaj u centar sfere
    kamera.ParallelProjectionOn()  # Ovo isključuje perspektivu (ortografski prikaz)

    # Inicijalizuj interakciju
    render_window.Render()
    
    # Parametri za animaciju
    trenutni_frejm = [0]  # Koristimo listu da bismo mogli da menjamo u callback-ovima
    ukupno_frejmova = len(fluid_aktori)
    vreme_izmedju_frejmova = int(1000 / fps)  # milisekunde
    
    # Timer callback za animaciju
    def timer_callback(caller, event):
        # Sakrij trenutni fluid aktor
        fluid_aktori[trenutni_frejm[0]].VisibilityOff()
        
        # Pređi na sledeći frejm
        trenutni_frejm[0] = (trenutni_frejm[0] + 1) % ukupno_frejmova
        
        # Prikaži novi fluid aktor
        fluid_aktori[trenutni_frejm[0]].VisibilityOn()
        
        # Ispiši informaciju o trenutnom frejmu
        print(f"Frejm: {trenutni_frejm[0] + 1}/{ukupno_frejmova} | Sfera na (3,2,3) r=0.35", end='\r')
        
        # Osveži prikaz
        render_window.Render()
    
    # Keyboard callback za kontrolu
    def keypress_callback(caller, event):
        key = caller.GetKeySym()
        
        if key == "space":  # Pauziraj/Nastavi
            if hasattr(interactor, 'GetTimerDuration') and interactor.GetTimerDuration() > 0:
                interactor.DestroyTimer()
                print("\nAnimacija pauzirana")
            else:
                interactor.CreateRepeatingTimer(vreme_izmedju_frejmova)
                print("\nAnimacija nastavljena")
        
        elif key == "Right":  # Sledeći frejm
            # Sakrij trenutni fluid aktor
            fluid_aktori[trenutni_frejm[0]].VisibilityOff()
            
            # Pređi na sledeći frejm
            trenutni_frejm[0] = (trenutni_frejm[0] + 1) % ukupno_frejmova
            
            # Prikaži novi fluid aktor
            fluid_aktori[trenutni_frejm[0]].VisibilityOn()
            render_window.Render()
            print(f"\nRučno: frejm {trenutni_frejm[0] + 1}/{ukupno_frejmova}")
        
        elif key == "Left":  # Prethodni frejm
            # Sakrij trenutni fluid aktor
            fluid_aktori[trenutni_frejm[0]].VisibilityOff()
            
            # Pređi na prethodni frejm
            trenutni_frejm[0] = (trenutni_frejm[0] - 1) % ukupno_frejmova
            
            # Prikaži novi fluid aktor
            fluid_aktori[trenutni_frejm[0]].VisibilityOn()
            render_window.Render()
            print(f"\nRučno: frejm {trenutni_frejm[0] + 1}/{ukupno_frejmova}")
        
        elif key == "r" or key == "R":  # Resetuj kameru
            renderer.ResetCamera()
            # Ponovo fokusiraj na sferu
            kamera = renderer.GetActiveCamera()
            kamera.SetFocalPoint(3, 2, 3)
            render_window.Render()
            print("\nKamera resetovana")
        
        elif key == "s" or key == "S":  # Prikaži/sakrij sferu
            if sfera_actor.GetVisibility():
                sfera_actor.VisibilityOff()
                print("\nSfera sakrivena")
            else:
                sfera_actor.VisibilityOn()
                print("\nSfera prikazana")
            render_window.Render()
        
        elif key == "q" or key == "Q":  # Izlaz
            interactor.ExitCallback()
    
    # Poveži callback-e
    interactor.AddObserver("TimerEvent", timer_callback)
    interactor.AddObserver("KeyPressEvent", keypress_callback)
    
    # Pokreni timer za animaciju
    interactor.CreateRepeatingTimer(vreme_izmedju_frejmova)
    
    print(f"\nAnimacija pokrenuta sa {fps} frejmova u sekundi")
    print("=" * 50)
    print("KONTROLE:")
    print("  SPACE - pauziraj/nastavi animaciju")
    print("  → - sledeći frejm (ručno)")
    print("  ← - prethodni frejm (ručno)")
    print("  R - resetuj kameru")
    print("  S - prikaži/sakrij sferu")
    print("  Q - izlaz")
    print("=" * 50)
    print(f"KOLIDER: Crvena sfera na (3, 2, 3), prečnik 0.7")
    print(f"FLUID: {len(fluid_aktori)} frejmova, plava boja (providnost 0.7)")
    print("-" * 50)
    
    # Pokreni interakciju
    interactor.Start()

def proveri_koliziju_sa_sferom(tacka, centar_sfere=(3, 2, 3), poluprecnik=0.35):
    """
    Proverava da li je tačka unutar sfere (kolidera).
    Korisno za debug - možete dodati u kod ako želite da proveravate koliziju.
    """
    x, y, z = tacka
    cx, cy, cz = centar_sfere
    
    rastojanje = np.sqrt((x - cx)**2 + (y - cy)**2 + (z - cz)**2)
    
    if rastojanje <= poluprecnik:
        return True  # Tačka je unutar kolidera
    else:
        return False  # Tačka je van kolidera

if __name__ == "__main__":
    # ZAMENITE OVU PUTANJU SA VAŠOM!
    putanja_do_foldera = r"./bin/results_obj"
    
    # Proveri da li folder postoji
    if not os.path.exists(putanja_do_foldera):
        print(f"Folder ne postoji: {putanja_do_foldera}")
        print("Molim vas, postavite tačnu putanju do vaših .obj fajlova")
        print("\nPrimer:")
        print('putanja_do_foldera = r"D:/Racun/Dimitrije/FluidFlow_NoTest/bin/results_obj"')
        
        # Pokušaj da nađe .obj fajlove u trenutnom folderu
        trenutni_folder = os.path.dirname(os.path.abspath(__file__))
        print(f"\nTražim .obj fajlove u trenutnom folderu: {trenutni_folder}")
        
        obj_fajlovi = glob.glob(os.path.join(trenutni_folder, "*.obj"))
        if obj_fajlovi:
            print(f"Pronađeno {len(obj_fajlovi)} .obj fajlova u trenutnom folderu!")
            putanja_do_foldera = trenutni_folder
            napravi_obj_animaciju_sa_koliderom(putanja_do_foldera, fps=10)
        else:
            print("Nema .obj fajlova ni u trenutnom folderu.")
    else:
        # Pokreni animaciju sa sferom
        napravi_obj_animaciju_sa_koliderom(putanja_do_foldera, fps=10)