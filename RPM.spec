
# This spec file can be used to build an RPM package for installation.
# **NOTE**
#     Version, Release, and tagname information should be updated for the
#     particular release to build an RPM for.

# ----------------------------------------------------------------------------
# Disable the creation of debug-info RPM package
# that contains stripped binary debug symbols
%global _enable_debug_package 0
%global debug_package %{nil}
# Disable stripped binary generation
%define __os_install_post %{nil}

# ----------------------------------------------------------------------------
Name:		espa-land-surface-temperature
Version:	0.0.1
Release:	1%{?dist}
Summary:	ESPA Land Surface Temperature Software

Group:		ESPA
License:	Nasa Open Source Agreement
URL:		https://github.com/USGS-EROS/espa-land-surface-temperature.git

BuildRoot:	%(mktemp -ud %{_tmppath}/%{name}-%{version}-%{release}-XXXXXX)
BuildArch:	x86_64
Packager:	USGS EROS LSRD

BuildRequires:	espa-common
Requires:	espa-common >= 1.4.0


# ----------------------------------------------------------------------------
%description
Provides science application executables for generating land surface temperature products.  This is a C and Python implementation which is statically built.


# ----------------------------------------------------------------------------
# Specify the repository tag/branch to clone and build from
%define tagname dev_v0.0.1
# Specify the name of the directory to clone into
%define clonedname %{name}-%{tagname}


# ----------------------------------------------------------------------------
%prep
# We don't need to perform anything here


# ----------------------------------------------------------------------------
%build

# Start with a clean clone of the repo
rm -rf %{clonedname}
git clone --depth 1 --branch %{tagname} %{url} %{clonedname}
# Build the applications
cd %{clonedname}
make BUILD_TYPE=-static


# ----------------------------------------------------------------------------
%install
# Start with a clean installation location
rm -rf %{buildroot}
# Install the applications for a specific path
cd %{clonedname}
make install PREFIX=%{buildroot}/usr/local

# ----------------------------------------------------------------------------
%clean
# Cleanup our cloned repository
rm -rf %{clonedname}
# Cleanup our installation location
rm -rf %{buildroot}


# ----------------------------------------------------------------------------
%files
%defattr(-,root,root,-)
# All sub-directories are automatically included
/usr/local/bin/*
/usr/local/%{name}


# ----------------------------------------------------------------------------
%changelog
* Wed Sep 03 2015 Ronald D Dilley <rdilley@usgs.gov>
- Build for Oct 2015 release
- Initial implementation
