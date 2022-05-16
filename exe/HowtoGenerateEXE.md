gui.py is the original file for gui dev. 
gui_exe.py is for creating .exe.
Follow steps from https://stackoverflow.com/a/60953781

- create resource_path() function
- use resource_path() for file opening:
  - for example, PhotoImage(
      file = resource_path("assets\\button_1.png"))
- run 'pyi-makespec --onefile gui_exe.py' to create the spec file for the compile and build process:
- in the generated gui_exe.spec file:
  - at top: added_files = [
    ("assets", "assets")
    ]
  - Then, change the line of datas=[], to datas=added_files,
  - Change hiddenimports to hiddenimports=["sklearn.metrics.pairwise", "sklearn.cluster", "sklearn.utils._typedefs", "sklearn.neighbors._partition_nodes"]
- go to 'build/exe' folder, run  pyinstaller ../gui_exe.spec
- find the generated .exe file in exe/dist. The exe file can be ran anywhere!
