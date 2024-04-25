function out = isAdmin
out = System.Security.Principal.WindowsPrincipal(...
        System.Security.Principal.WindowsIdentity.GetCurrent()).IsInRole(...
        System.Security.Principal.WindowsBuiltInRole.Administrator);
end