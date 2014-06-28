Module cfml_python
USE CFML_IO_Formats
Use CFML_Crystal_Metrics,           only: Crystal_Cell_Type
Use CFML_Crystallographic_Symmetry, only: Space_Group_Type
Use CFML_Atom_TypeDef,              only: atom_list_type
CONTAINS
SUBROUTINE ReadXtal_Structure_File(filenam, Cell, SpG, A, Mode, Iphase, Job_Info, file_list, CFrame)
!---- Arguments ----!
       character(len=*), intent(in) :: filenam
       Type (Crystal_Cell_Type), intent(out) :: Cell
       Type (Space_Group_Type), intent(out) :: SpG
       Type (atom_list_type), intent(out) :: A
       Character(len=*), optional, intent(in) :: Mode
       Integer, optional, intent(in) :: Iphase
       Type(Job_Info_type), optional, intent(out) :: Job_Info
       Type(file_list_type), optional, intent(in out) :: file_list
       Character(len=*), optional, intent(in) :: CFrame
       call Readn_Set_Xtal_Structure(filenam, Cell, SpG, A, Mode, Iphase, Job_Info, file_list, CFrame)
END SUBROUTINE ReadXtal_Structure_File
End Module cfml_python
