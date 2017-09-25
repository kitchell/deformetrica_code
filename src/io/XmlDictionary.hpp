/***************************************************************************************
*                                                                                      *
*                                     Deformetrica                                     *
*                                                                                      *
*    Copyright Inria and the University of Utah.  All rights reserved. This file is    *
*    distributed under the terms of the Inria Non-Commercial License Agreement.        *
*                                                                                      *
*                                                                                      *
****************************************************************************************/

#ifndef DEFORMETRICA_XMLDICTIONARY_H
#define DEFORMETRICA_XMLDICTIONARY_H

#include "itkObject.h"
#include "itkDOMReader.h"
#include "itkDOMNodeXMLReader.h"
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>
#include <typeinfo>
#include <functional>
#include <algorithm>
#include <memory>
#include <map>

namespace def {
namespace io {

namespace filters {

struct filter {
  ~filter() {}
  virtual std::string operator()(const std::string &v) = 0;
};

struct lower_case : public filter {
  virtual std::string operator()(const std::string &v) {
    return itksys::SystemTools::LowerCase(v);
  }
};

}

namespace range {
enum value {
  negative_include_zero,
  negative_exclude_zero,
  positive_include_zero,
  positive_exclude_zero
};
}

class XmlProgrammable {
 public:
  typedef std::function<void(const std::string &val)> reader;

  XmlProgrammable() {
    assign_ = [](const std::string &val) {};
    range_ = one_of_ = assign_;
  }

  XmlProgrammable& required(bool r) {
    required_ = r;
    return *this;
  }

  const bool is_required() const {
    return required_;
  }

  template <typename T>
  XmlProgrammable& add_reader(T f) {
    readers_.push_back(f);
    return *this;
  }
//enable_if<is_pointer_t<T>::value, void>::type
//  template <typename T, typename F, typename O,typename std::is_same<decltype(O::operator*), decltype(O::operator->)>::value>
  template <typename T, typename O, typename F>
  XmlProgrammable& assign_to(O o, F f) {
    assign_ = [o, f](const std::string &val) {
      std::stringstream ss;
      T data;
      ss << val;
      ss >> data;
      if (o == nullptr)
        throw std::runtime_error("assing_to() methond has null object");
      (o->*f)(data);
    };

    return *this;
  }

  template <typename F>
  XmlProgrammable& filter_with(F f) {
    filter_ = std::make_shared<F>();
    return *this;
  }

  template <class T, typename... Args>
  XmlProgrammable& one_of(T first, Args... args) {
    auto candidates = std::vector<T>({first, args...});

    if (typeid(T()) == typeid(std::string()))
      for (auto&it :candidates)
        it = itksys::SystemTools::LowerCase(it);

    one_of_ = [candidates] (const std::string &val) {
      std::stringstream ss;
      T data;
      ss << val;
      ss >> data;

      if (typeid(T()) == typeid(std::string()))
        data = itksys::SystemTools::LowerCase(data);

      if (std::find(candidates.begin(), candidates.end(), data) != candidates.end())
        return;

      ss.clear();
      ss << "The input value [" << val << "] is not valid..." << std::endl;
      ss << "Allowed data are : ";
      std::copy(candidates.begin(), candidates.end(), std::ostream_iterator<T>(ss," "));

      throw std::domain_error(ss.str());
    };

    return *this;
  }

  template <typename T>
  XmlProgrammable& range(T min, T max) {
    range_ = [min,max](const std::string& val) {
      std::stringstream ss;
      T data;
      ss << val;
      ss >> data;

      if (data < min || data > max)
        throw std::domain_error(std::string("Range not allowed for value [" + val + "]"));
    };

    return *this;
  }

  template <typename T>
  XmlProgrammable& range(range::value mask) {
    range_ = [mask](const std::string& val) {
      std::stringstream ss;
      T data;
      ss << val;
      ss >> data;

      switch (mask) {
        case range::value::negative_exclude_zero:
          if (data >= 0)
            throw std::domain_error(std::string("Positive and zero number is not allowed!"));
          break;
        case range::value::negative_include_zero:
          if (data > 0)
            throw std::domain_error(std::string("Positive number is not allowed!"));
          break;
        case range::value::positive_exclude_zero:
          if (data <= 0)
            throw std::domain_error(std::string("Negative and zero number is not allowed!"));
          break;
        case range::value::positive_include_zero:
          if (data < 0)
            throw std::domain_error(std::string("Negative number is not allowed!"));
          break;
      }
    };

    return *this;
  }

  XmlProgrammable& range(reader&& r) {
    range_ = r;

    return *this;
  }

  void run(const std::string &val) {
    std::string v = val;
    if (filter_)
      v = filter_->operator()(val);

    one_of_(v);
    range_(v);

    for (auto reader : readers_)
      reader(v);

    assign_(v);
  }

 protected:
  std::shared_ptr<filters::filter> filter_;
  reader assign_;
  reader one_of_;
  reader range_;
  std::vector<reader> readers_;
  bool required_;
};


class XmlDictionary : public XmlProgrammable {
 private:
  typedef std::map<std::string, XmlDictionary> XmlMap;
  typedef std::map<std::string, std::shared_ptr<XmlMap>> TagMap;
  XmlMap map_;
  TagMap tag_;

 public:
  XmlDictionary() : XmlProgrammable() { required_ = false; }

  XmlDictionary& operator[](std::string key) {
    if (map_.find(key) == map_.end()) {
      XmlDictionary t;
      map_[key] = std::move(t);
    }

    return map_[key];
  }

  typedef XmlMap::iterator iterator;
  typedef XmlMap::const_iterator const_iterator;

  iterator begin() { return  map_.begin(); }
  iterator end() { return  map_.end(); }

  const_iterator cbegin() { return map_.cbegin(); }
  const_iterator cend() { return map_.cend(); }

  iterator begin(const std::string& tag) { return  tag_[tag]->begin(); }
  iterator end(const std::string& tag) { return  tag_[tag]->end(); }

  const_iterator cbegin(const std::string& tag) { return tag_[tag]->cbegin(); }
  const_iterator cend(const std::string& tag) { return tag_[tag]->cend(); }

  XmlDictionary& add_tag(const std::string& tag_id, const std::string& tag_value) {
    if ((tag_.find(tag_id) == tag_.end())) {
      XmlMap m;
      XmlDictionary t;
      m.emplace(tag_value, std::move(t));
      tag_[tag_id] = std::make_shared<XmlMap>(m);
    } else if (tag_[tag_id]->find(tag_value) == tag_[tag_id]->end()) {
      XmlDictionary t;
      tag_[tag_id]->emplace(tag_value, t);
    }

    return tag_[tag_id]->at(tag_value);
  }

  const bool load_tag(const std::string& tag_id, const std::string& tag_value, XmlDictionary& xml) {
    if (tag_.find(tag_id) == tag_.end()) {
      return false;
    }

    if (tag_[tag_id]->find(tag_value) == tag_[tag_id]->end()) {
      return false;
    }

    xml = tag_[tag_id]->at(tag_value);
    return true;
  }

  const bool has_tags() const {
    return (tag_.size() > 0);
  }

  std::vector<std::string> tags() {
    std::vector<std::string> tag;
    for (auto it : tag_)
      tag.push_back(it.first);

    return tag;
  }

  decltype(map_.size()) size() const { return map_.size(); }

  static void xml_reader(itk::DOMNode::Pointer dom, XmlDictionary &main_dictionary, std::string first_tag=std::string()) {
    static const auto cmp = [](const std::string& v1, const std::string& v2) {
      return (itksys::SystemTools::LowerCase(v1) == itksys::SystemTools::LowerCase(v2));
    };

    if (!first_tag.empty()) {
      if (!cmp(dom->GetName(), first_tag))
        throw std::runtime_error(std::string("Expected <" + first_tag + "> as first xml tag"));
    }

    itk::DOMNode::ChildrenListType all_children;
    dom->GetAllChildren(all_children);

    for (auto it : main_dictionary) {
      auto &key_dictionary = it.first;
      auto &sub_dictionary = it.second;
      itk::DOMNode::ChildrenListType children;

      //resolve case-insensitive
      for (auto ch : all_children) {
        if (!cmp(ch->GetName(),key_dictionary)) continue;

        children.push_back(ch);
      }

      if(!children.size()) {
        if (sub_dictionary.is_required())
          throw std::runtime_error(std::string("Tag <" + key_dictionary + "> is required and is missing"));

        continue;
      }

      if (sub_dictionary.tags().size()) {
        for (auto id_tag : sub_dictionary.tags()) {
          for (auto xml_tag = sub_dictionary.begin(id_tag); xml_tag != sub_dictionary.end(id_tag); ++xml_tag) {
            bool find_tag = false;
            itk::DOMNode::Pointer dom_child;
            for (auto it : children) {
              if (cmp(it->GetAttribute(id_tag), xml_tag->first)) {
                dom_child = it;//->Clone();
                find_tag = true;
                break;
              }
            }

            if (/*xml_tag->second.is_required() && */!find_tag)
              throw std::runtime_error(std::string("Tag <" + key_dictionary + " " + id_tag + "=\"" + xml_tag->first + "\"> is missing"));

            if (find_tag) {
              xml_reader(dom_child, xml_tag->second);
            }
          }
        }
      }

      if(children.size() == 1) {
        auto dom_child = children[0];
        if (sub_dictionary.size()) {
          xml_reader(dom_child, sub_dictionary);
        } else {
          if (dom_child->GetTextChild()) {
            try {
              auto str = dom_child->GetTextChild()->GetText();
              sub_dictionary.run(str);
            } catch (std::exception& ex) {
              std::stringstream ss;
              ss << std::endl << "Error working on tag <" << dom_child->GetName() << ">" << std::endl;
              ss << ex.what();
              throw std::runtime_error(ss.str());
            }
          }
        }
      }
    }
  }

};


}
}

#endif //DEFORMETRICA_XMLDICTIONARY_H
