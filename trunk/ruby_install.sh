#!/bin/bash

# First thisngs first, try to find what version of ruby we have
INSTALL_RUBY=0
RUBY_WEBSITE="http://ftp.ruby-lang.org/pub/ruby"
RUBY_V=`ruby -v`
# capture the exit value
EV=$?

if test $EV == 127; then
	INSTALL_RUBY=1
	RUBY_V="ruby-1.9.3-p0"
	RUBY_MAJMIN="1.9"
else
	# check the version here, if it's too old, we'll need to get a newer
	# version
	RUBY_MAJMIN=`echo $RUBY_V | sed "s/ruby \([0-9]\.[0-9]\).*/\1/g"`
	RUBY_V=`echo $RUBY_V | sed "s/ruby \([0-9]\.[0-9]\.[0-9]\)p\([0-9]*\).*/ruby-\1-p\2/g"`
fi

if test $INSTALL_RUBY -eq 0; then
	# download and install ruby
	echo "Installing ruby v1.9.3p0..."
	curl "$RUBY_WEBSITE/$RUBY_MAJMIN/$RUBY_V.tar.gz" > $RUBY_V.tar.gz
	echo "Downloaded"
	tar -xzf $RUBY_V.tar.gz
	cd  $RUBY_V
	# I'm assuming that you have access to the necessary development
	# tools to compile ruby from source
	./configure
	YAML_INST = `make 2>&1 | grep "libyaml is missing"`
	echo "YAML $YAML_INST"
	exit 1
	if test !-z $YAML_INST; then
		echo "LibYAML not found, installing..."
		# We need to install libyaml first
		curl http://pyyaml.org/download/libyaml/yaml-0.1.4.tar.gz > yaml-0.14.tar.gz
		tar -xzf yaml-0.1.4.tar.gz
		cd yaml-0.1.4
		./configure 
		make
		make install
		cd ..
		make clean
		make
	fi
	make install
	cd ..
fi

INSTALL_GEMS=0
INSTALL_RMAGICK=0
GEM_LIST=`gem list`
EV=$?

if test $EV == 127; then
	INSTALL_GEMS=1
fi

if test -z "`echo $GEM_LIST | grep rmagick`"; then
	INSTALL_RMAGICK=1
fi

if test $INSTALL_GEMS -eq 1; then
	echo "Installing RubyGems v1.8.12..."
	# download and install rubygems for the appropriate version here
	curl http://production.cf.rubygems.org/rubygems/rubygems-1.8.12.tgz > rubygems-1.8.12.tgz
	tar -xzf rubygems-1.8.12.tgz
	cd rubygems-1.8.12
	ruby setup.rb
	cd ..	
fi

if test $INSTALL_RMAGICK -eq 1; then
	echo "Installing RMagick..."
	# Use rubygems to install rmagick
	gem install rmagick
fi

#OK, you should be good to go now!
