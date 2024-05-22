import os, re, stat, sys, shutil

from PyQt5 import QtWidgets, QtCore

qtWidgets_modules = dir(QtWidgets)
qtCore_modules = dir(QtCore)

def fixFile(line_list):
	line_generator = iter(line_list)
	changed = False
	num = -1

	for each in line_generator:
		num += 1
		if 'from PyQt4 import' in each.strip():
			line_list[num] = 'from PyQt5 import QtWidgets, QtGui, QtCore\n'
			changed = True
			continue

		result = re.search('QtGui\.(\w+)', each)

		if result and result.group(1) in qtWidgets_modules:
			line_list[num] = each.replace('QtGui', 'QtWidgets')
			changed = True

		if result and result.group(1) in qtCore_modules:
			line_list[num] = each.replace('QtGui', 'QtCore')
			changed = True

	if changed:
		return ''.join(line_list)

def convert(main_path):
	for dirs, folders, files in os.walk(main_path):
		for each in [x for x in files if x.endswith('.py')]:
			file_path = os.path.join(dirs, each).replace(os.sep, '/')
			os.chmod(file_path, stat.S_IWRITE)

			with open(file_path, 'r+') as out_file:
				line_list = out_file.readlines()

			fixed_lines = fixFile(line_list)

			if fixed_lines:
				## Copy and rename file into a folder name old_qt4 so source is not lost
				dir_name, file_name = os.path.split(file_path)
				new_location = os.path.join(dir_name, 'old_qt4', file_name)
				shutil.copy(file_path, new_location)
				
				## Rewrite file with changes from Qt4 to Qt5
				with open(file_path, 'w+') as out_file:
					out_file.write(fixed_lines)

if __name__ == "__main__":
	## Get directory path from user
	folder_path = input('Where is the directory that needs to be converted:')
	
	## Check if path is a directory
	if os.path.isdir(folder_path) :
		convert(main_path=folder_path)
		print('All ".py" files in %s have been converted to Qt5')
	else:
		print('File path: %s does not exists'%(folder_path))