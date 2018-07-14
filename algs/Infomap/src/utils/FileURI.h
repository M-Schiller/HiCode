/**********************************************************************************

 Infomap software package for multi-level network clustering

 Copyright (c) 2013 Daniel Edler, Martin Rosvall

 For more information, see <http://www.mapequation.org>

 This file is part of Infomap software package.

 Infomap software package is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 Infomap software package is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with Infomap software package.  If not, see <http://www.gnu.org/licenses/>.

**********************************************************************************/
#pragma once

#ifndef FILEURI_H_
#define FILEURI_H_

#include <iostream>
#include <string>

/**
 * Filename class to simplify handling of parts of a filename.
 * If a path is path/to/file.ext, the member methods of this class give these parts:
 * getFilename -> "path/to/file.ext"
 * getDirectory -> "path/to/"
 * getName -> "file"
 * getExtension -> "ext"
 */
class FileURI
{
public:
  FileURI(); // Good to be able to have a Filename member in a class without require initialization in initialization list
  FileURI(const char* filename, bool requireExtension = false);
  FileURI(const std::string& filename, bool requireExtension = false);
  FileURI(const FileURI& other);
  FileURI & operator= (const FileURI& other);

  const std::string& getFilename() const
  {
    return m_filename;
  }

  /**
   * Includes last '/' if non-empty.
   */
  const std::string& getDirectory() const
  {
    return m_directory;
  }

  const std::string& getName() const
  {
    return m_name;
  }

  const std::string& getExtension() const
  {
    return m_extension;
  }

  std::string getParts() const
  {
    std::string out("['" +
      getDirectory() + "' + '" +
      getName() + "' + '" +
      getExtension() + "']");
    return out;
  }

  friend std::ostream& operator<<(std::ostream& out, const FileURI& file)
  {
    return out << file.getFilename();
  }

private:
  void analyzeFilename();
  std::string getErrorMessage();

  std::string m_filename;
  bool m_requireExtension;

  std::string m_directory;
  std::string m_name;
  std::string m_extension;
};

#endif /* FILEURI_H_ */
